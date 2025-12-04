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

        public double From(double x, double y, double z)
        {
            double px = x, py = y, pz = z;
            double minDistSq = double.MaxValue;

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
                    double distSq = dx * dx + dy * dy + dz * dz;
                    if (distSq < minDistSq)
                        minDistSq = distSq;
                }
                else if (t >= 1.0)
                {
                    // Closest to v2
                    double dx = px - x2;
                    double dy = py - y2;
                    double dz = pz - z2;
                    double distSq = dx * dx + dy * dy + dz * dz;
                    if (distSq < minDistSq)
                        minDistSq = distSq;
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

                    double distSq = dx * dx + dy * dy + dz * dz;
                    if (distSq < minDistSq)
                        minDistSq = distSq;
                }
            }

            return Math.Sqrt(minDistSq);
        }
    }
}