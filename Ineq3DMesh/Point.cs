using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;

namespace MeshData
{
    public class Point : IComparable
    {
        public static bool EnableUnsafeMove = false;

        private double x = 0;
        private double y = 0;
        private double z = 0;
        private double u = 0;
        private BitArray boundary;

        private HashSet<Tetrahedron> tetrahedrons = new HashSet<Tetrahedron>();

        private Point() { }

        public Point(double x, double y, double z)
        {
            this.x = x;
            this.y = y;
            this.z = z;
        }

        public Point(double x, double y, double z, int ineqCount)
        {
            this.x = x;
            this.y = y;
            this.z = z;
            boundary = new BitArray(ineqCount);
        }
        
        public double X
        {
            get { return x; }
            set { x = value; }
        }

        public double Y
        {
            get { return y; }
            set { y = value; }
        }

        public double Z
        {
            get { return z; }
            set { z = value; }
        }

        public double U
        {
            get { return u; }
            set { u = value; }
        }

        public BitArray Boundary
        {
            get { return boundary; }
        }

        public int BoundaryCount
        {
            get { return  Boundary.Cast<bool>().Count(b => b); }
        }

        public int BoundaryFirstIndex
        {
            get
            {
                return Boundary.Cast<bool>().Select((b, i) => new { IsSet = b, index = i }).First(bf => bf.IsSet).index;
            }
        }

        public int BoundarySecondIndex
        {
            get
            {
                return Boundary.Cast<bool>().Select((b, i) => new { IsSet = b, index = i }).Where(bf => bf.IsSet).Skip(1).First().index;
            }
        }

        public int BoundaryThirdIndex
        {
            get
            {
                return Boundary.Cast<bool>().Select((b, i) => new { IsSet = b, index = i }).Where(bf => bf.IsSet).Skip(2).First().index;
            }
        }

        public bool HasAllBoundary(Point p)
        {
            return Boundary.Cast<bool>().Select((b, i) => new { IsSet = b, index = i }).Where(bf => bf.IsSet).All(bf => p.Boundary[bf.index]);
        }

        public int[] HasNotBoundary(Point p)
        {
            var res = Boundary.Cast<bool>().Select((b, i) => new { IsSet = b, index = i }).Where(bf => bf.IsSet && !p.Boundary[bf.index]).Select(bf => bf.index).ToArray();
            return res;
        }


        public bool HasAnyBoundary(Point p)
        {
            return Boundary.Cast<bool>().Select((b, i) => new { IsSet = b, index = i }).Where(bf => bf.IsSet).Any(bf => p.Boundary[bf.index]);
        }

        public HashSet<Tetrahedron> Tetrahedrons
        {
            get { return tetrahedrons; }
        }

        public IEnumerable<Point> Points
        {
            get 
            { 
                return tetrahedrons
                    .SelectMany(t => t.Points)
                    .Where(p => p != this)
                    .Distinct(); 
            }
        }

        public bool IsOnBoundary
        {
            get
            {
                return Tetrahedrons
                    .SelectMany(t => t.Edges())
                    .Where(p => p.P1 != this && p.P2 != this)
                    .GroupBy(p => p)
                    .Any(g => g.Count() == 1);
            }
        }


        public override string ToString()
        {
            return string.Format("[{0}; {1}; {2}]", x.ToString("n2"), y.ToString("n2"), z.ToString("n2"));
        }

        public double Distance(Point p)
        {
            return Math.Sqrt((p.x - x) * (p.x - x) + (p.y - y) * (p.y - y) + (p.z - z) * (p.z - z));
        }

        public static Point operator +(Point p1, Point p2)
        {
            return new Point(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z, p1.boundary.Length);
        }

        public static Point operator -(Point p1, Point p2)
        {
            return new Point(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z, p1.boundary.Length);
        }

        public static Point operator /(Point p, double f)
        {
            return new Point(p.x / f, p.y / f, p.z / f, p.boundary.Length);
        }

        public bool MoveTo(Point p1, bool safe)
        {
            return MoveTo(p1.X, p1.Y, p1.Z, safe);
        }
        
        public bool MoveTo(double x, double y, double z, bool safe)
        {
            if (EnableUnsafeMove || !safe)
            {
                this.x = x; this.y = y; this.z = z;

                return true;
            }
            else
            {
                double origX = this.x;
                double origY = this.y;
                double origZ = this.z;

                this.x = x; this.y = y; this.z = z;

                int p = 0;
                while(p < 20 && Tetrahedrons.Any(t => !t.CheckVolume()))
                {
                    this.x = origX + 0.9 * (this.x - origX);
                    this.y = origY + 0.9 * (this.y - origY);
                    this.z = origZ + 0.9 * (this.z - origZ);

                    p++;
                }

                if (p == 20)
                {
                    this.x = origX;
                    this.y = origY;
                    this.z = origZ;
                    return false;
                }
                else
                    return true;
                
            }
        }

        public bool CanMoveTo(Point p1)
        {
            return CanMoveTo(p1.X, p1.Y, p1.Z);
        }

        public bool CanMoveTo(double x, double y, double z)
        {
            double origX = this.x;
            double origY = this.y;
            double origZ = this.z;

            bool res = true;

            this.x = x; this.y = y; this.z = z;

            if (Tetrahedrons.Any(t => !t.CheckVolume()))
            {
                res = false;
            }

            this.x = origX;
            this.y = origY;
            this.z = origZ;

            return res;
        }

        public bool CanMoveTo(Point p1, IEnumerable<Tetrahedron> exceptTetras, double tolerance = 0.001d)
        {
            double origX = this.x;
            double origY = this.y;
            double origZ = this.z;

            bool res = true;

            this.x = p1.X; this.y = p1.Y; this.z = p1.Z;

            if (Tetrahedrons.Except(exceptTetras).Any(t => !t.CheckVolume(tolerance)))
            {
                res = false;
            }

            this.x = origX;
            this.y = origY;
            this.z = origZ;

            return res;
        }



        //!!!!! GetHashCode method does not guarantee unique return values for different objects, replace with instance-unique value !!!!!
        public static bool operator >(Point p1, Point p2)
        {
            return  p1.GetHashCode() > p2.GetHashCode();
        }

        public static bool operator <(Point p1, Point p2)
        {
            return !(p1.GetHashCode() > p2.GetHashCode());
        }

        public int CompareTo(object obj)
        {
            return this.GetHashCode().CompareTo(obj.GetHashCode());
        }
    }
}
