using System;
using System.Threading;
using System.Threading.Tasks;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace MeshData
{
    public class Mesh
    {
        protected HashSet<Point> points = new HashSet<Point>();
        protected HashSet<Tetrahedron> tetrahedrons = new HashSet<Tetrahedron>();

        protected int boundaryCount = 0;
        protected int domainCount = 0;

        public void Clear()
        {
            points = new HashSet<Point>();
            tetrahedrons = new HashSet<Tetrahedron>(); 
        }
        
        public HashSet<Point> Points
        {
            get { return points; }
        }

        public HashSet<Tetrahedron> Tetrahedrons
        {
            get { return tetrahedrons; }
        }

        public virtual Point AddPoint(double x, double y, double z)
        {
            Point p = new Point(x, y, z);
            points.Add(p);
            return p;
        }

        public Point AddPoint(Point p)
        {
            points.Add(p);
            return p;
        }

        public void DeletePoint(Point p)
        {
            foreach(Tetrahedron t in p.Tetrahedrons.ToArray())
            {
                DeleteTetrahedron(t);
            }
            points.Remove(p);
        }

        public Tetrahedron AddTetrahedron(Point p0, Point p1, Point p2, Point p3, int volumeSign = -1)
        {
            Tetrahedron t = new Tetrahedron(p0, p1, p2, p3, boundaryCount, domainCount, false, volumeSign);
            tetrahedrons.Add(t);
            return t;
        }

        public void DeleteTetrahedron(Tetrahedron t)
        {
            foreach (Point p in t.Points)
            {
                p.Tetrahedrons.Remove(t);
            } 
            tetrahedrons.Remove(t);
        }

        public ParallelQuery<Edge> Edges
        {
            get
            {
                return Points.AsParallel()//.WithExecutionMode(ParallelExecutionMode.ForceParallelism)
                    .SelectMany(
                        p1 => p1.Tetrahedrons.SelectMany(t => t.Points.Where(p2 => p2 > p1)).Distinct(),
                        (p1, p2) => new Edge(p1, p2)
                        );
            }
        }

        public ParallelQuery<Edge> EdgesFromPoints(IEnumerable<Point> fromPoints)
        {
            return fromPoints.AsParallel()//.WithExecutionMode(ParallelExecutionMode.ForceParallelism)
                .SelectMany(
                    p1 => p1.Tetrahedrons.SelectMany(t => t.Points.Where(p2 => p2 > p1)).Distinct(),
                    (p1, p2) => new Edge(p1, p2)
                    );
                //.Distinct();
        }

        public ParallelQuery<Edge> EdgesForBoundary(int ineqNumber)
        {

            return Points.AsParallel()//.WithExecutionMode(ParallelExecutionMode.ForceParallelism)
                    .SelectMany(
                        p1 => p1.Tetrahedrons.Where(t => t.Boundary[ineqNumber]).SelectMany(t => t.Points.Where(p2 => p2 > p1)).Distinct(),
                        (p1, p2) => new Edge(p1, p2)
                        );
        }


        public IEnumerable<Triangle> Triangles
        {
            get
            {
                return Tetrahedrons.SelectMany(t => t.Triangles());
            }
        }

        public IEnumerable<Triangle> TrianglesForTetrahedrons(IEnumerable<Tetrahedron> tetrahedrons)
        {
            return tetrahedrons.SelectMany(t => t.Triangles());
        }


        public IEnumerable<Point> Boundary
        {
            get
            {
                return Points.Where(p => p.IsOnBoundary);
            }
        }



        public IEnumerable<Point> LonelyPoints
        {
            get
            {
                return Points.AsParallel()
                    .Where(p => !p.Tetrahedrons.Any());
            }
        }

        public void DeleteLonelyPoints()
        {
            foreach (Point p in LonelyPoints.ToArray())
            {
                DeletePoint(p);
            } 
        }
    }

    public struct Edge : IEnumerable<Point>
    {
        Point p1;
        Point p2;

        public double SqrLength
        {
            get
            {
                return (p1.X - p2.X) * (p1.X - p2.X) + (p1.Y - p2.Y) * (p1.Y - p2.Y) + (p1.Z - p2.Z) * (p1.Z - p2.Z);
            }
        }
        
        public Edge(Point p1, Point p2)
        {
            if (p1 < p2)
            {
                this.p1 = p1; this.p2 = p2;
            }
            else
            {
                this.p2 = p1; this.p1 = p2;
            }
        }

        public Point P1
        {
            get { return p1; }
        }

        public Point P2
        {
            get { return p2; }
        }

        public bool Valid
        {
            get
            {
                Point p2 = this.p2;
                return p1.Tetrahedrons.Any(t => t.Points.Contains(p2));
            }
        }

        public IEnumerator<Point> GetEnumerator()
        {
            yield return p1;
            yield return p2;
        }

        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator()
        {
            yield return p1;
            yield return p2;
        }
    }


    public struct Triangle : IEnumerable<Point>
    {
        Point p1;
        Point p2;
        Point p3;


        public Triangle(Point p1, Point p2, Point p3/*, bool sort = true*/)
        {
            Point[] points = new Point[3];

            points[0] = p1;
            points[1] = p2;
            points[2] = p3;

            /*if (sort)
            {
                Array.Sort(points);
            }*/

            this.p1 = points[0];
            this.p2 = points[1];
            this.p3 = points[2];   
        }

        public Point P1
        {
            get { return p1; }
        }

        public Point P2
        {
            get { return p2; }
        }

        public Point P3
        {
            get { return p3; }
        }

        public int BoundaryCount
        {
            get
            {
                return 
                    P1.Boundary.Cast<bool>().Select((b, i) => new { b = b, i = i }).Where(bi => bi.b).Select(bi => bi.i)
                .Intersect
                (
                    P2.Boundary.Cast<bool>().Select((b, i) => new { b = b, i = i }).Where(bi => bi.b).Select(bi => bi.i)
                )
                .Intersect
                (
                    P3.Boundary.Cast<bool>().Select((b, i) => new { b = b, i = i }).Where(bi => bi.b).Select(bi => bi.i)
                ).Count();
            }
        }

        public Nullable<int> CommonBoundaryFlag
        {
            get
            {
                return
                    P1.Boundary.Cast<bool>().Select((b, i) => new { b = b, i = i }).Where(bi => bi.b).Select(bi => bi.i)
                .Intersect
                (
                    P2.Boundary.Cast<bool>().Select((b, i) => new { b = b, i = i }).Where(bi => bi.b).Select(bi => bi.i)
                )
                .Intersect
                (
                    P3.Boundary.Cast<bool>().Select((b, i) => new { b = b, i = i }).Where(bi => bi.b).Select(bi => bi.i)
                ).Cast<int?>().FirstOrDefault();
            }
        }

        public bool Boundary
        {
            get
            {
                return P1.Tetrahedrons.Intersect(P2.Tetrahedrons).Intersect(P3.Tetrahedrons).Count() == 1;
            }
        }

        public IEnumerator<Point> GetEnumerator()
        {
            yield return p1;
            yield return p2;
            yield return p3;
        }

        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator()
        {
            yield return p1;
            yield return p2;
            yield return p3;
        }

        public IEnumerable<Edge> Edges()
        {
            yield return new Edge(p1, p2);
            yield return new Edge(p1, p3);
            yield return new Edge(p2, p3);
        }

        public bool Valid
        {
            get
            {
                Point p2 = this.p2;
                Point p3 = this.p3;
                return p1.Tetrahedrons.Intersect(p2.Tetrahedrons).Intersect(p3.Tetrahedrons).Any();
            }
        }

        public double[]  Normal()
        {
            Point u = p2 - p1;
            Point v = p3 - p1;

            double n1 = u.Y * v.Z - u.Z * v.Y;
            double n2 = -(u.X * v.Z - u.Z * v.X);
            double n3 = u.X * v.Y - u.Y * v.X;
            double n = Math.Sqrt(n1 * n1 + n2 * n2 + n3 * n3);

            if (n != 0)
            {
                return new double[] { n1 / n, n2 / n, n3 / n };
            }
            else
            {
                return new double[] {0, 0, 0};
            }
        }

        public double Quality
        {
            get
            {
                double l = Math.Sqrt(Edges().Max(e => e.SqrLength));
                double vMax = Math.Sqrt(3) * l * l / 4;

                Point u = p2 - p1;
                Point v = p3 - p1;

                double n1 = u.Y * v.Z - u.Z * v.Y;
                double n2 = -(u.X * v.Z - u.Z * v.X);
                double n3 = u.X * v.Y - u.Y * v.X;
                double a = Math.Sqrt(n1 * n1 + n2 * n2 + n3 * n3) / 2;

                return a / vMax;
            }
        }


        public override bool Equals(object obj)
        {
            if (!(obj is Triangle other))
                return false;

            Point[] points = new Point[3];
            Point[] points1 = new Point[3];

            points[0] = p1;
            points[1] = p2;
            points[2] = p3;

            points1[0] = ((Triangle)obj).P1;
            points1[1] = ((Triangle)obj).P2;
            points1[2] = ((Triangle)obj).P3;

            Array.Sort(points);
            Array.Sort(points1);

            return points[0].Equals(points1[0]) && points[1].Equals(points1[1]) & points[2].Equals(points1[2]);
        }

        public static bool operator ==(Triangle left, Triangle right)
        {
            return left.Equals(right);
        }

        public static bool operator !=(Triangle left, Triangle right)
        {
            return !(left == right);
        }

        public override int GetHashCode()
        {
            unchecked
            {
                Point[] points = new Point[3];

                points[0] = p1;
                points[1] = p2;
                points[2] = p3;

                Array.Sort(points);

                int hash = 17; 
                hash = hash * 23 + points[0].GetHashCode(); 
                hash = hash * 23 + points[1].GetHashCode(); 
                hash = hash * 23 + points[2].GetHashCode();
                return hash;
            }
        }

        public class IntersectionResult
        {
            public string Type;   // "inside", "edge", "vertex", or "none"
            public Point Point;   // intersection point (if any)
            public Edge? Edge;   // optional: "P1P2", "P2P3", "P3P1"
        }

        public IntersectionResult IntersectLine(Point A, Point B, double eps = 1e-10)
        {
            var AB = B - A;
            var N = Point.Cross(P2 - P1, P3 - P1);  // plane normal
            double denom = Point.Dot(N, AB);

            // Line and plane parallel?
            if (Math.Abs(denom) < eps)
                return new IntersectionResult { Type = "none", Point = null };

            // Compute intersection parameter t
            double t = Point.Dot(N, P1 - A) / denom;
            var I = A + t * AB; // intersection point

            // Now check if I is inside triangle (using barycentric coords)
            var v0 = P2 - P1;
            var v1 = P3 - P1;
            var v2 = I - P1;

            double d00 = Point.Dot(v0, v0);
            double d01 = Point.Dot(v0, v1);
            double d11 = Point.Dot(v1, v1);
            double d20 = Point.Dot(v2, v0);
            double d21 = Point.Dot(v2, v1);

            double denomBary = d00 * d11 - d01 * d01;
            if (Math.Abs(denomBary) < eps)
                return new IntersectionResult { Type = "none", Point = null }; // degenerate triangle

            double v = (d11 * d20 - d01 * d21) / denomBary;
            double w = (d00 * d21 - d01 * d20) / denomBary;
            double u = 1.0 - v - w;

            // Check if inside or on edge/vertex
            if (u < -eps || v < -eps || w < -eps)
                return new IntersectionResult { Type = "none", Point = null };

            // Classify location
            string type;
            Edge? edge = null;

            bool onU = Math.Abs(u) < eps;
            bool onV = Math.Abs(v) < eps;
            bool onW = Math.Abs(w) < eps;

            int onCount = (onU ? 1 : 0) + (onV ? 1 : 0) + (onW ? 1 : 0);

            if (onCount == 3)
                type = "vertex"; // Degenerate triangle actually
            else if (onCount == 2)
            {
                type = "vertex";
                if (onU && onV) I = P3;
                else if (onV && onW) I = P1;
                else if (onW && onU) I = P2;
            }
            else if (onCount == 1)
            {
                type = "edge";
                if (onU) edge = new Edge(P2, P3);// "P2P3";
                else if (onV) edge = new Edge(P3, P1); //"P3P1";
                else if (onW) edge = new Edge(P1, P2);//"P1P2";
            }
            else
                type = "inside";

            return new IntersectionResult { Type = type, Point = I, Edge = edge };
        }
    }


    public static class MeshTools
    {
        public static Point Average(this IEnumerable<Point> points)
        {
            double x = 0; double y = 0; double z = 0;
            int count = 0;
            foreach (Point p in points)
            {
                count++;
                x += p.X;
                y += p.Y;
                z += p.Z;
            }
            if (count > 0)
                return new Point(x / count, y / count, z / count);
            else
                return null;
        }
    }
}
