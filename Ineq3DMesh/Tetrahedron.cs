using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;

namespace MeshData
{
    public class Tetrahedron
    {
        private Point[] points = new Point[4];
        private BitArray isin;
        private BitArray isOnBoundary;
        private BitArray boundary;

        private double origVolume = 0;

        private Tetrahedron() { }

        internal Tetrahedron(Point p0, Point p1, Point p2, Point p3, int boundaryCount, int domainCount, bool standAlone = false, int volumeSign = -1)
        {
            points[0] = p0;
            points[1] = p1;
            points[2] = p2;
            points[3] = p3;

            if (!standAlone)
            {
                //Array.Sort(points);
                
                p0.Tetrahedrons.Add(this);
                p1.Tetrahedrons.Add(this);
                p2.Tetrahedrons.Add(this);
                p3.Tetrahedrons.Add(this);

                isin = new BitArray(domainCount/*2 * p0.Boundary.Length - 1*/);
                isOnBoundary = new BitArray(domainCount/*2 * p0.Boundary.Length - 1*/);
                boundary = new BitArray(boundaryCount/*p0.Boundary.Length*/);
            }

            OrigVolume = Volume;

            if (Math.Sign(OrigVolume) != volumeSign)
            {
                var tmp = points[0];
                points[0] = points[1];
                points[1] = tmp;
                OrigVolume = -OrigVolume;
            }
        }

        public Point[] Points
        {
            get { return points; }
        }

        public Point P0
        {
            get
            {
                return points[0];
            }
        }

        public Point P1
        {
            get
            {
                return points[1];
            }
        }

        public Point P2
        {
            get
            {
                return points[2];
            }
        }

        public Point P3
        {
            get
            {
                return points[3];
            }
        }


        public BitArray IsIn
        {
            get { return isin; }
        }

        public BitArray IsOnBoundary
        {
            get { return isOnBoundary; }
        }

        public bool IsInDomain(int domaiNumber)
        {
            return IsIn[domaiNumber] && !IsOnBoundary[domaiNumber];
        }

        public bool IsOutDomain(int domaiNumber)
        {
            return !IsIn[domaiNumber] && !IsOnBoundary[domaiNumber];
        }

        public bool IsOnBoundaryDomain(int domaiNumber)
        {
            return !IsIn[domaiNumber] && IsOnBoundary[domaiNumber];
        }



        public BitArray Boundary
        {
            get { return boundary; }
        }

        public int BoundaryCount
        {
            get { return Boundary.Cast<bool>().Count(b => b); }
        }

        public IEnumerable<Edge> Edges()
        {
            yield return new Edge(P0, P1);
            yield return new Edge(P0, P2);
            yield return new Edge(P0, P3);


            yield return new Edge(P1, P2);
            yield return new Edge(P1, P3);

            yield return new Edge(P2, P3);
        }

        public IEnumerable<Triangle> Triangles()
        {
            yield return new Triangle(P0, P1, P2/*, false*/); //3
            yield return new Triangle(P1, P0, P3/*, false*/); //2
            yield return new Triangle(P0, P2, P3/*, false*/); //1
            yield return new Triangle(P2, P1, P3/*, false*/); //0
        }


        public void ReplacePoint(Point oldPoint, Point newPoint)
        {
            for (int i = 0; i < 4; i++)
            {
                if (points[i] == oldPoint)
                {
                    points[i] = newPoint;
                    oldPoint.Tetrahedrons.Remove(this);
                    newPoint.Tetrahedrons.Add(this);
                }
            }
        }

        public double OrigVolume
        {
            get
            {
                return origVolume;
            }
            private set
            {
                origVolume = value;
           }
        }
        
        public double Volume
        {
            get
            {
                Point v1 = P1 - P0;
                Point v2 = P2 - P0;
                Point v3 = P3 - P0;

                return Numeric.Determinant(v1.X, v1.Y, v1.Z, v2.X, v2.Y, v2.Z, v3.X, v3.Y, v3.Z);
            }
        }

        public double Quality
        {
            get
            {
                double l = Math.Sqrt(Edges().Max(e => e.SqrLength));
                double vMax = Math.Sqrt(2) * l * l * l / 12;
                double v = Math.Abs(Volume) / 6;

                return v / vMax;
            }
        }

        public bool CheckVolume(double tolerance = 0.01d)
        {
            double volume = Volume;
            return OrigVolume * volume > 0 && Math.Abs(volume) > tolerance * Math.Abs(OrigVolume); 
        }

        public int CommonBoundaryCount
        {
            get
            {
                return
                    P0.Boundary.Cast<bool>().Select((b, i) => new { b = b, i = i }).Where(bi => bi.b).Select(bi => bi.i)
                .Intersect
                (
                    P1.Boundary.Cast<bool>().Select((b, i) => new { b = b, i = i }).Where(bi => bi.b).Select(bi => bi.i)
                )
                .Intersect
                (
                    P2.Boundary.Cast<bool>().Select((b, i) => new { b = b, i = i }).Where(bi => bi.b).Select(bi => bi.i)
                )
                .Intersect
                (
                    P3.Boundary.Cast<bool>().Select((b, i) => new { b = b, i = i }).Where(bi => bi.b).Select(bi => bi.i)
                )                
                .Count();
            }
        }
    }
}
