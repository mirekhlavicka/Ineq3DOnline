using MeshData;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Web;

namespace Ineq3DOnline
{
    public class WireFrameDistance
    {
        private class Vertex
        {
            public Double X;
            public Double Y;
            public Double Z;
        }

        private class Edge
        {
            public Vertex v1;
            public Vertex v2;
        }

        //private List<Vertex> vertices = new List<Vertex>();
        private List<Edge> edges = new List<Edge>();

        public void AddEdge(double x1, double y1, double z1, double x2, double y2, double z2)
        {
            var v1 = new Vertex { X = x1, Y = y1, Z = z1 };
            var v2 = new Vertex { X = x2, Y = y2, Z = z2 };
            edges.Add(new Edge {v1 = v1, v2 = v2});
        }

        public double From(double x, double y, double z, double eps = 0)
        {
            double px = x, py = y, pz = z;
            double minDist = 1000000;//double.MaxValue;

            foreach (var e in edges)
            {
                // Edge endpoints
                double x1 = e.v1.X, y1 = e.v1.Y, z1 = e.v1.Z;
                double x2 = e.v2.X, y2 = e.v2.Y, z2 = e.v2.Z;

                // Vector v = v2 - v1
                double vx = x2 - x1;
                double vy = y2 - y1;
                double vz = z2 - z1;

                // Vector w = P - v1
                double wx = px - x1;
                double wy = py - y1;
                double wz = pz - z1;

                double vv = vx * vx + vy * vy + vz * vz; // |v|^2

                double t;
                double dist;
                if (vv < 1e-20)
                {
                    // v1 and v2 are the same point → treat as point distance
                    t = 0.0;
                }
                else
                {
                    // Projection parameter t = dot(w, v) / dot(v, v)
                    t = (wx * vx + wy * vy + wz * vz) / vv;
                }

                if (t <= 0.0)
                {
                    // Closest to v1
                    double dx = wx;
                    double dy = wy;
                    double dz = wz;
                    dist = Math.Sqrt(dx * dx + dy * dy + dz * dz);
                }
                else if (t >= 1.0)
                {
                    // Closest to v2
                    double dx = px - x2;
                    double dy = py - y2;
                    double dz = pz - z2;
                    dist = Math.Sqrt(dx * dx + dy * dy + dz * dz);
                }
                else
                {
                    // Projection point is between v1 and v2
                    double projx = x1 + t * vx;
                    double projy = y1 + t * vy;
                    double projz = z1 + t * vz;

                    double dx = px - projx;
                    double dy = py - projy;
                    double dz = pz - projz;

                    dist = Math.Sqrt(dx * dx + dy * dy + dz * dz);
                }

                if (eps == 0)
                {
                    if (dist < minDist)
                        minDist = dist;
                }
                else
                {

                    minDist = SmoothMin(minDist, dist, eps);
                }
            }

            return minDist;
        }

        private double SmoothMin(double a, double b, double eps)
        {
            return (a + b - MeshData.IneqTree.IneqNode.SmoothAbs(a - b, eps)) / 2.0d;
        }

        public IneqTree ToPrism(double r, int count, double f = 1.0d, double[] v = null, double sa = 0)
        { 
            var res = new IneqTree();

            foreach (var e in edges)
            {
                double x1 = e.v1.X, y1 = e.v1.Y, z1 = e.v1.Z;
                double x2 = e.v2.X, y2 = e.v2.Y, z2 = e.v2.Z;
                double xc = (x1 + x2) / 2;
                double yc = (y1 + y2) / 2;
                double zc = (z1 + z2) / 2;

                x1 = xc + f * (x1 - xc); y1 = yc + f * (y1 - yc); z1 = zc + f * (z1 - zc);
                x2 = xc + f * (x2 - xc); y2 = yc + f * (y2 - yc); z2 = zc + f * (z2 - zc);

                res = res | IneqLib.Prism(x1, y1, z1, x2, y2, z2, r, count, v, sa);
            }

            return res;
        }

        /*public static double SmoothAbs(double x, double eps)
        {

            return Math.Sqrt(x * x + eps * eps);
        }*/
    }
}