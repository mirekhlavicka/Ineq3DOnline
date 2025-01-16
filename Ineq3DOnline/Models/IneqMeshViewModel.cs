using MeshData;
using System;
using System.Collections.Generic;
using System.ComponentModel.DataAnnotations;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Web;

namespace Ineq3DOnline.Models
{
    public class IneqMeshViewModel
    {

        private static string[] sampleFormulas =
        {
@"z < (sin(6*x)+sin(6*y))/4",
//------------------------------------------------------------
@"(
	x^2 + y^2 + z^2 < 0.75^2
) ||
(
	x^2 + y^2 + (z-0.75)^2 < 0.25^2
)",
//------------------------------------------------------------
@"(
    x^2 + y^2 + z^2 < 1.25 ^ 2
) &&
x^2+y^2 > 0.5^2 &&
x^2+z^2 > 0.5^2 &&
z^2+y^2 > 0.5^2",
//------------------------------------------------------------
@"(
    x^2 + y^2 + z^2 < 1 &&
    x + y < z^2
) ||
(
    x^2 + z^2 < 0.25
)",
//------------------------------------------------------------
@"(
	x^2+y^2+z^2<1 &&
	z<0+(cos(5*x)+cos(5*y))/3
) ||
x^2+y^2+(z-.8)^2<0.25",
//------------------------------------------------------------
@"abs(x)<0.25 ||
0.5*x^2+y^2+z^2<0.5",
//------------------------------------------------------------
@"8*x^2+8*y^2+z^2<1 ||
8*x^2+y^2+8*z^2<1 ||
x^2+8*y^2+8*z^2<1",
//------------------------------------------------------------
@"16*x^2+16*y^2+z^2<1 ||
90*x^6+y^2+90*z^6<1 ||
x^2+8*y^2+8*z^2<1",
//------------------------------------------------------------
@"(x-0.25*cos(4*z))^2+(y-0.25*sin(4*z))^2<0.85 &&
(x-0.25*cos(4*z))^2+(y-0.25*sin(4*z))^2>0.3"

        };

        private static string[] sampleUFuncs = {
            "sin(6 * x - t) + sin(6 * y - t) + sin(6 * z - t)",
            "x * y * (z - sin(t)^2)",
            "sin(6 * x - 0.1*t^2) + sin(6 * y - t) + sin(6 * z - t)",
            "(x-sin(t))^2+y^2+z^2",
            "(x-sin(t))^2+cos(t)*y^2+z^2",
            "x * y ^2 * (z - sin(t)) +0.4*y*(x-cos(2*t))^3+z",
            "sin(5*((x-sin(t))^2+cos(t)*y^2+z^2))",
            "sin(4*x-4*sin(t))+sin(4*y)+sin(4*z)",
            "sin(t)*(x^2+y^2+z^2)+(1-sin(t))*(x*y*z)",
            "sin(t)*(z-x^2+y^2)+(1-sin(t))*(x-y^2-z^2)"
        } ;


        public int SampleIndex { get; set; }
        public int SampleUFuncIndex { get; set; }

        public static IneqMeshViewModel DefaultModel(int sample = -1, int sampleU = -1)
        {
            if (sample == -1)
            {
                sample = new Random().Next(sampleFormulas.Length);
            }
            else
            {
                sample = sample % sampleFormulas.Length;
            }

            if (sampleU == -1)
            {
                sampleU = new Random().Next(sampleUFuncs.Length);
            }
            else
            {
                sampleU = sampleU % sampleUFuncs.Length;
            }

            return new IneqMeshViewModel
            {
                Formula = sampleFormulas[sample],
                SampleIndex = sample,
                SampleUFuncIndex = sampleU,
                MaxDivisionCount = 12,
                X0 = -1,
                Y0 = -1,
                Z0 = -1,
                X1 = 1,
                Y1 = 1,
                Z1 = 1,
                UFunc = sampleUFuncs[sampleU]
            };
        }

        [Required(ErrorMessage = "The Inequalities are required.")]
        public string Formula { get; set; }

        [Required(ErrorMessage = "The Division count is required.")]
        [Range(5, 25, ErrorMessage = "Value for {0} must be between {1} and {2}.")]
        [Display(Name = "Division count")]
        public int MaxDivisionCount { get; set; }

        public double X0 { get; set; }
        public double Y0 { get; set; }
        public double Z0 { get; set; }
        public double X1 { get; set; }
        public double Y1 { get; set; }
        public double Z1 { get; set; }

        [Display(Name = "Improve mesh quality")]
        public bool Quality { get; set; }

        [Display(Name = "Improve approximation by curvature")]
        public bool CurvatureQuality { get; set; }


        //[Required(ErrorMessage = "The Function is required.")]
        [Display(Name = "u(x,y,z,t) = ")]
        public string UFunc { get; set; }

        public IneqMesh IneqMesh { get; set; }

        public void SetIneqMesh()
        {
            IneqMesh = new IneqMesh
            {
                X0 = X0,
                Y0 = Y0,
                Z0 = Z0,
                X1 = X1,
                Y1 = Y1,
                Z1 = Z1,
                D = Math.Max(Math.Max(X1 - X0, Y1 - Y0), Z1 - Z0) / MaxDivisionCount,
                Boxed = true,
                IneqTree = IneqTreeParser.FromFormula(Formula)
            };
        }

        public void CreateMesh(bool createPLY = true) 
        {
            if (IneqMesh == null)
            {
                return;
            }

            IneqMesh.Create();

            if (Quality)
            {
                CheckQuality(IneqMesh);
                CheckBoundaryQuality(IneqMesh);
            }

            if (CurvatureQuality)
            {
                CheckCurvatureQuality(IneqMesh);
            }
            
            IneqMesh.DeleteLonelyPoints();

            if (createPLY)
            {
                PLY = GetPLY(IneqMesh);
            }
        }

        double minQuality = 0.25d;
        public void CheckQuality(IneqMesh ineqMesh)
        {
            if (ineqMesh == null)
                return;

            int c = 0;
            Dictionary<int, int> counts = new Dictionary<int, int>();
            do
            {
                c = ineqMesh.CheckQuality(minQuality, false);

                if (!counts.Keys.Contains(c))
                    counts[c] = 1;
                else
                    counts[c]++;


            }
            while (c != 0 && counts[c] < 3);

            if (c != 0)
            {
                c = ineqMesh.CheckQuality(minQuality, true);
            }

            ineqMesh.DeleteLonelyPoints();

            ineqMesh.Jiggle(3);

            return;
        }

        double bminQuality = 0.25d;
        public void CheckBoundaryQuality(IneqMesh ineqMesh)
        {
            if (ineqMesh == null)
                return;

            int c = 0;
            Dictionary<int, int> counts = new Dictionary<int, int>();
            do
            {
                c = ineqMesh.CheckBoundaryQuality(minQuality, false);

                if (!counts.Keys.Contains(c))
                    counts[c] = 1;
                else
                    counts[c]++;

            }
            while (c != 0 && counts[c] < 3);

            //if (c != 0)
            //{
            //    c = ineqMesh.CheckBoundaryQuality(minQuality, true);
            //}

            ineqMesh.DeleteLonelyPoints();

            ineqMesh.Jiggle(3);

            return;
        }

        public void CheckCurvatureQuality(IneqMesh ineqMesh)
        {
            List<Triangle> refList = new List<Triangle>();

            var edges = ineqMesh.Edges
                .Where(e => e.P1.Boundary.Cast<bool>().Select((b, i) => new { b = b, i = i }).Where(bi => bi.b && e.P2.Boundary[bi.i]).Count() == 2)
                .Select(e =>
                {
                    var cf = e.P1.Boundary.Cast<bool>().Select((b, i) => new { b = b, i = i }).Where(bi => bi.b && e.P2.Boundary[bi.i]).ToArray();

                    return new
                    {
                        bf1 = cf[0].i,
                        bf2 = cf[1].i,
                        p = e.Average(),
                        e = e,
                        length = e.P1.Distance(e.P2)
                    };
                }).ToArray();

            foreach (var ee in edges)
            {
                if (!ee.e.Valid)
                {
                    continue;
                }
                Point origp = new Point(ee.p.X, ee.p.Y, ee.p.Z);

                ineqMesh.ProjectToEdge(ee.p, ee.bf1, ee.bf2, false);

                double dist = origp.Distance(ee.p);

                if (dist >= ee.length / 25.0d)
                {
                    var e = ee.e;

                    var trians = e
                        .P1.Tetrahedrons.Intersect(e.P2.Tetrahedrons)
                        .SelectMany(tt => tt.Triangles())
                        .Where(tr => tr.Contains(e.P1) && tr.Contains(e.P2))
                        .GroupBy(tr => tr)
                        .Where(gr => gr.Count() == 1)
                        .Select(gr => gr.Single());

                    refList.AddRange(trians);

                }
            }

            var centerPoints = ineqMesh.Tetrahedrons.SelectMany(t => t.Triangles()
                            .Where(tr => tr.BoundaryCount == 1 && tr.P1.Tetrahedrons.Intersect(tr.P2.Tetrahedrons).Intersect(tr.P3.Tetrahedrons).Count() == 1))
                            .Select(tr => new
                            {
                                bf = tr.CommonBoundaryFlag.Value,
                                p = tr.Average(),
                                tr = tr,
                                maxLength = Math.Max(Math.Max(tr.P1.Distance(tr.P2), tr.P1.Distance(tr.P3)), tr.P2.Distance(tr.P3))
                            });

            foreach (var cp in centerPoints)
            {
                Point origp = new Point(cp.p.X, cp.p.Y, cp.p.Z);

                ineqMesh.ProjectToSurface(cp.p, 100, cp.bf, false);

                double dist = origp.Distance(cp.p);

                if (dist >= cp.maxLength / 25.0d)
                {
                    refList.Add(cp.tr);
                }
            }

            ineqMesh.RefineBoundaryTriangles(refList);

            ineqMesh.Jiggle(3);
        }

        public string PLY { get; set; }

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
        private string[] colors = {
            "255 0 0",
            "0 255 0",
            "0 0 255",
            "255 255 0",
            "0 255 255",
            "255 0 255",
            "255 255 255",
            "255 255 128",
            "128 255 255",
            "255 128 255"
        };

        public string GetPLY(IneqMesh ineqMesh)
        {
            StringBuilder sbPoints = new StringBuilder();
            StringBuilder sbTriangles = new StringBuilder();
            int pointsCount = 0;

            var boundaryTriangles = ineqMesh.Tetrahedrons.SelectMany(t => t.Triangles().Where(tr => /*tr.BoundaryCount == 1 &&*/ tr.P1.Tetrahedrons.Intersect(tr.P2.Tetrahedrons).Intersect(tr.P3.Tetrahedrons).Count() == 1))
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

        public byte[] ConvertPLYToSTL()
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