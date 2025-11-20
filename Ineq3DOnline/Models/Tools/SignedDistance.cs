using MeshData;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Web;
using System.Xml.Linq;

namespace Ineq3DOnline
{
	public class SignedDistance
	{
        private TriangleMeshInterop.TriangleMeshWrapper meshWrapper = null;
        List<double[]> v = null;
        List<int[]> t = null;

        List<double[]> V 
        {
            get
            {
                return v;
            }
        }

        List<int[]> T
        {
            get
            {
                return t;
            }
        }

        public SignedDistance(string filePath, double d = 2, double x0 = 0, double y0 = 0, double z0 = 0)
		{

            filePath = System.IO.Path.Combine(HttpContext.Current.Server.MapPath("~/MeshSamples"), filePath);

            STLLoader.LoadSTL(filePath, out v, out t);

            /*List<int[]> bad = new List<int[]>();

            foreach (var tr in t)
            {
                if (!CheckTriangle(v[tr[0]], v[tr[1]], v[tr[2]]))
                {
                    bad.Add(tr);
                }
            }

            foreach (var tr in bad)
            {
                t.Remove(tr);
            }*/

            t.RemoveAll(tr => !CheckTriangle(v[tr[0]], v[tr[1]], v[tr[2]]));

            if (d > 0)
            {
                double minX = double.MaxValue, minY = double.MaxValue, minZ = double.MaxValue;
                double maxX = -double.MaxValue, maxY = -double.MaxValue, maxZ = -double.MaxValue;

                foreach (var p in v)
                {
                    minX = Math.Min(minX, p[0]);
                    minY = Math.Min(minY, p[1]);
                    minZ = Math.Min(minZ, p[2]);
                    maxX = Math.Max(maxX, p[0]);
                    maxY = Math.Max(maxY, p[1]);
                    maxZ = Math.Max(maxZ, p[2]);
                }

                double s = Math.Min(Math.Min(d / (maxX - minX), d / (maxY - minY)), d / (maxZ - minZ));
                double x00 = (maxX + minX) / 2.0d;
                double y00 = (maxY + minY) / 2.0d;
                double z00 = (maxZ + minZ) / 2.0d;

                foreach (var p in v)
                {
                    p[0] = x0 + s * (p[0] - x00);
                    p[1] = y0 + s * (p[1] - y00);
                    p[2] = z0 + s * (p[2] - z00);
                }
            }

            meshWrapper = new TriangleMeshInterop.TriangleMeshWrapper(v, t);
        }

        private static bool CheckTriangle(double[] v1, double[] v2, double[] v3)
        {
            double[] u = { v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2] };
            double[] v = { v3[0] - v1[0], v3[1] - v1[1], v3[2] - v1[2] };

            double nx = u[1] * v[2] - u[2] * v[1];
            double ny = u[2] * v[0] - u[0] * v[2];
            double nz = u[0] * v[1] - u[1] * v[0];

            return Math.Sqrt(nx * nx + ny * ny + nz * nz) > 1e-6; //0

            /*if (Math.Sqrt(nx * nx + ny * ny + nz * nz) > 1e-6)
            {
                return true;
            }
            else
            {
                return false;
            }*/            
        }

        public double From(double x, double y, double z)
        {
            double[] result = meshWrapper.SignedDistance(x, y, z);
            return result[0];
        }

        public bool Project(Point p)
        {
            double[] result = meshWrapper.SignedDistance(p.X, p.Y, p.Z);

            p.X = result[1];
            p.Y = result[2];
            p.Z = result[3];

            return true;
        }

        public List<Tuple<int, int, double>> FindSharpEdges(double thresholdDegrees)
        {
            double threshold = thresholdDegrees * Math.PI / 180.0;

            // edge → triangles
            var edgeTris = new Dictionary<Tuple<int, int>, List<int>>();

            for (int ti = 0; ti < t.Count; ti++)
            {
                int[] tri = t[ti];
                AddEdge(edgeTris, tri[0], tri[1], ti);
                AddEdge(edgeTris, tri[1], tri[2], ti);
                AddEdge(edgeTris, tri[2], tri[0], ti);
            }

            var sharpEdges = new List<Tuple<int, int, double>>();

            foreach (var kv in edgeTris)
            {
                var edge = kv.Key;
                List<int> tris = kv.Value;

                if (tris.Count != 2)
                    continue; // skip boundary edges

                int t1 = tris[0];
                int t2 = tris[1];

                double angle = DihedralAngle(v, t[t1], t[t2], edge.Item1, edge.Item2);

                double deviation = Math.Abs(angle - Math.PI);

                if (deviation > threshold)
                {
                    sharpEdges.Add(Tuple.Create(edge.Item1, edge.Item2, angle));
                }
            }

            return sharpEdges;
        }

        private static void AddEdge(
            Dictionary<Tuple<int, int>, List<int>> dict,
            int a, int b, int triIndex)
        {
            if (a > b) { int tmp = a; a = b; b = tmp; }

            var key = Tuple.Create(a, b);

            List<int> list;
            if (!dict.TryGetValue(key, out list))
            {
                list = new List<int>();
                dict[key] = list;
            }
            list.Add(triIndex);
        }

        private static double DihedralAngle(
            List<double[]> v,
            int[] tri1,
            int[] tri2,
            int a, int b)
        {
            var A = v[a];
            var B = v[b];

            int c1 = ThirdVertex(tri1, a, b);
            int c2 = ThirdVertex(tri2, a, b);

            var C1 = v[c1];
            var C2 = v[c2];

            var n1 = Normalize(Cross(Sub(B, A), Sub(C1, A)));
            var n2 = Normalize(Cross(Sub(B, A), Sub(C2, A)));

            double dot = Dot(n1, n2);
            if (dot < -1) dot = -1;
            if (dot > 1) dot = 1;

            return Math.Acos(dot);
        }

        private static int ThirdVertex(int[] tri, int a, int b)
        {
            if (tri[0] != a && tri[0] != b) return tri[0];
            if (tri[1] != a && tri[1] != b) return tri[1];
            return tri[2];
        }

        // -------- vector helpers --------

        private static double[] Sub(double[] a, double[] b)
        {
            return new double[] { a[0] - b[0], a[1] - b[1], a[2] - b[2] };
        }

        private static double Dot(double[] a, double[] b)
        {
            return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        }

        private static double[] Cross(double[] a, double[] b)
        {
            return new double[]
            {
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
            };
        }

        private static double[] Normalize(double[] v)
        {
            double len = Math.Sqrt(Dot(v, v));
            return new double[] { v[0] / len, v[1] / len, v[2] / len };
        }

        public void ForceSharpEdges(IneqMesh ineqMesh, IneqTree ineq, double thresholdDegrees)
        {
            foreach (var e in FindSharpEdges(thresholdDegrees))
            {
                var p1 = new Point(v[e.Item1][0], v[e.Item1][1], v[e.Item1][2]);
                var p2 = new Point(v[e.Item2][0], v[e.Item2][1], v[e.Item2][2]);

                ineqMesh.ForceEdge(p1, p2, ineq);
            }
        }
    }
}