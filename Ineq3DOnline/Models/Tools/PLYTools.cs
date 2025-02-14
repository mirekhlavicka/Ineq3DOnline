using MeshData;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Web;

namespace Ineq3DOnline
{
    public static class PLYTools
    {
        const string plyFormat =
@"ply
format ascii 1.0
element vertex {0}
property float32 x
property float32 y
property float32 z
property uchar red
property uchar green
property uchar blue
property uchar alpha
element face {1}
property list uint8 int32 vertex_indices
end_header
{2}
{3}
";
        private static string[] colors = {
            "255 128 128",
            "128 255 128",
            "128 128 255",

            "128 255 255",
            "255 128 255",
            "255 255 128",

            "128 255 0",
            "0 128 255",
            "255 0 128",

            "0 255 255",
            "255 0 255",
            "255 255 0",

            "255 0 0",
            "0 255 0",
            "0 0 255",

            "255 255 255"
        };

        public static string GetPLY(IneqMesh ineqMesh)
        {
            StringBuilder sbPoints = new StringBuilder();
            StringBuilder sbTriangles = new StringBuilder();
            int pointsCount = 0;

            var boundaryTriangles = ineqMesh.Tetrahedrons.SelectMany(t => t.Triangles().Where(tr => /*tr.BoundaryCount == 1 &&*/ tr.Boundary))
                //.Where(t => t.Quality > 0.25d)
                .ToArray();

            for (int ineqNumber = 0; ineqNumber < ineqMesh.IneqTree.ExpressionList.Count + 3; ineqNumber++)
            {
                var triangles = boundaryTriangles.Where(t => t.CommonBoundaryFlag == ineqNumber).ToArray();
                var points = triangles.SelectMany(t => t).Distinct().Select((p, i) => new { Point = p, index = i }).ToDictionary(pi => pi.Point, pi => pi.index);

                foreach (var p in points.OrderBy(p => p.Value))
                {
                    sbPoints.Append(String.Format(System.Globalization.CultureInfo.InvariantCulture, "{0} {1} {2} {3} 255 {4}", p.Key.X, p.Key.Y, p.Key.Z, colors[ineqNumber % colors.Length], System.Environment.NewLine));
                }

                foreach (var t in triangles)
                {
                    sbTriangles.Append(String.Format(System.Globalization.CultureInfo.InvariantCulture, "3 {0} {1} {2} {3}", points[t.P1] + pointsCount, points[t.P2] + pointsCount, points[t.P3] + pointsCount, System.Environment.NewLine));
                }

                pointsCount += points.Count;
            }

            return String.Format(plyFormat, pointsCount, boundaryTriangles.Length, sbPoints, sbTriangles);
        }

        public static byte[] ConvertPLYToSTL(string PLY)
        {
            string plyContent = PLY;
            List<double[]> vertices = new List<double[]>();
            List<int[]> faces = new List<int[]>();
            bool isHeader = true;
            bool hasVertex = false;
            bool hasFace = false;
            int vertexCount = 0;
            int faceCount = 0;

            using (StringReader sr = new StringReader(plyContent))
            {
                string line;
                while ((line = sr.ReadLine()) != null)
                {
                    if (String.IsNullOrEmpty(line))
                    {
                        continue;
                    }

                    if (isHeader)
                    {
                        if (line.StartsWith("element vertex"))
                        {
                            vertexCount = int.Parse(line.Split(' ')[2]);
                            hasVertex = true;
                        }
                        else if (line.StartsWith("element face"))
                        {
                            faceCount = int.Parse(line.Split(' ')[2]);
                            hasFace = true;
                        }
                        else if (line.StartsWith("end_header"))
                        {
                            isHeader = false;
                        }
                    }
                    else
                    {
                        if (hasVertex && vertices.Count < vertexCount)
                        {
                            var vertex = Array.ConvertAll(line.Split(' ').Take(3).ToArray(), s => double.Parse(s, CultureInfo.InvariantCulture));
                            vertices.Add(vertex);
                        }
                        else if (hasFace && faces.Count < faceCount)
                        {
                            var face = Array.ConvertAll(line.Split(' ').Skip(1).Take(3).ToArray(), int.Parse);
                            faces.Add(face);
                        }
                    }
                }
            }

            return GetBinarySTL(vertices, faces);
        }

        private static byte[] GetBinarySTL(List<double[]> vertices, List<int[]> faces)
        {
            using (MemoryStream ms = new MemoryStream())
            using (BinaryWriter bw = new BinaryWriter(ms))
            {
                // 80-byte header
                bw.Write(new byte[80]);

                // Number of triangles
                bw.Write(faces.Count);

                foreach (var face in faces)
                {
                    var v1 = vertices[face[0]];
                    var v2 = vertices[face[1]];
                    var v3 = vertices[face[2]];

                    var normal = CalculateNormal(v1, v2, v3);
                    bw.Write((float)normal[0]);
                    bw.Write((float)normal[1]);
                    bw.Write((float)normal[2]);

                    foreach (var vertex in new[] { v1, v2, v3 })
                    {
                        bw.Write((float)vertex[0]);
                        bw.Write((float)vertex[1]);
                        bw.Write((float)vertex[2]);
                    }

                    bw.Write((ushort)0); // Attribute byte count
                }

                return ms.ToArray();
            }
        }

        private static double[] CalculateNormal(double[] v1, double[] v2, double[] v3)
        {
            double[] u = { v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2] };
            double[] v = { v3[0] - v1[0], v3[1] - v1[1], v3[2] - v1[2] };

            double nx = u[1] * v[2] - u[2] * v[1];
            double ny = u[2] * v[0] - u[0] * v[2];
            double nz = u[0] * v[1] - u[1] * v[0];

            double length = Math.Sqrt(nx * nx + ny * ny + nz * nz);
            return new double[] { nx / length, ny / length, nz / length };
        }
    }
}