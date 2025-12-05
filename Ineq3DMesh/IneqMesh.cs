using System;
using System.Threading;
using System.Threading.Tasks;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;

namespace MeshData
{
    using FuncXYZ = Func<double, double, double, double>;

    public class IneqMesh : Mesh
    {
        private const double maxWarpDist = 0.5d;// 0.48d 

        public double X0 { get; set; }
        public double Y0 { get; set; }
        public double Z0 { get; set; }
        public double X1 { get; set; }
        public double Y1 { get; set; }
        public double Z1 { get; set; }

        public double D { get; set; }

        public Action<double> OnProgress { get; set; }

        public Dictionary<Func<double, double, double, double>, Func<Point, bool>> ProjectToSurfaceSpec = new Dictionary<Func<double, double, double, double>, Func<Point, bool>>();

        private IneqTree ineqTreeBoxed;
        private IneqTree ineqTree;

        private int currentDomaiNumber = 0;

        public IneqTree IneqTree
        {
            get { return ineqTree; }
            set { ineqTree = value; }
        }

        public IneqTree IneqTreeBoxed
        {
            get { return ineqTreeBoxed; }
        }

        public Action PrepareBackgroundMesh { get; set; }

        public Action PrepareBackgroundMeshBeforeApriory { get; set; }

        public bool Boxed { get; set; }

        public void Create()
        {
            Point.EnableUnsafeMove = false;
            Clear();

            if (Boxed)
            {
                ineqTreeBoxed = ineqTree &
                    ((IneqTree)
                    ((x, y, z) => Math.Abs(x - (X0 + X1) / 2.0d) - 0.5d * (X1 - X0)) &
                    ((x, y, z) => Math.Abs(y - (Y0 + Y1) / 2.0d) - 0.5d * (Y1 - Y0)) &
                    ((x, y, z) => Math.Abs(z - (Z0 + Z1) / 2.0d) - 0.5d * (Z1 - Z0)));
            }
            else
                ineqTreeBoxed = ineqTree;

            boundaryCount = ineqTreeBoxed.ExpressionList.Count;
            domainCount = ineqTreeBoxed.Root.DomainCount;//2 * boundaryCount - 1;


            CreateBackgroundMesh();
            DeleteLonelyPoints();

            if (PrepareBackgroundMeshBeforeApriory != null)
            {
                PrepareBackgroundMeshBeforeApriory();
                //JiggleBackgroundMash(5);
            }

            ResolveMeshApriori(ineqTreeBoxed.Root, 0);
            SpreadTetrahedronsBoundaryFlags();

            /*if (PrepareBackgroundMesh != null)
            {
                SpreadTetrahedronsBoundaryFlags();
            }*/

            if (PrepareBackgroundMesh != null)
            {
                PrepareBackgroundMesh();
            }

            foreach (Tetrahedron t in Tetrahedrons.AsParallel().Where(t => t.BoundaryCount == 0 && !t.IsIn[0]).ToArray())
                DeleteTetrahedron(t);
            DeleteLonelyPoints();

            ResolveMesh(ineqTreeBoxed.Root, 0);
            foreach (Tetrahedron t in Tetrahedrons.AsParallel().Where(t => !t.IsIn[0]).ToArray())
                DeleteTetrahedron(t);
            DeleteLonelyPoints();

            //var c = Points.Count(p => p.BoundaryCount >= 3);
            //CheckTopology();
            //c = Points.Count(p => p.BoundaryCount >= 3);

            Jiggle(3);

            if (PrepareBackgroundMesh == null && PrepareBackgroundMeshBeforeApriory == null)
            {
                CheckQuality(0.25d, false);
                CheckQuality(0.25d, false);
                CheckTopology();
            }

            Point.EnableUnsafeMove = true;
            Jiggle(3);
        }

        private void ResolveMesh(IneqTree.IneqNode node, int domaiNumber)
        {
            if (domaiNumber == 0)
                currentDomaiNumber = 0;

            if (node.NodeType == IneqTree.NodeType.NodeExpression)
            {
                ResolveIneq(node.ExpressionIndex, domaiNumber);
            }
            else
            {
                int domaiNumberLeft = ++currentDomaiNumber;
                int domaiNumberRight = ++currentDomaiNumber;

                ResolveMesh(node.Left, domaiNumberLeft);
                ResolveMesh(node.Right, domaiNumberRight);

                Parallel.ForEach(Tetrahedrons, t =>
                {
                    if (node.NodeType == IneqTree.NodeType.NodeOr)
                    {
                        t.IsIn[domaiNumber] = t.IsIn[domaiNumberLeft] || t.IsIn[domaiNumberRight];
                    }
                    if (node.NodeType == IneqTree.NodeType.NodeAnd)
                    {
                        t.IsIn[domaiNumber] = t.IsIn[domaiNumberLeft] && t.IsIn[domaiNumberRight];
                    }
                });

                int[] expressionIndexes = node.ExpressionIndexes.ToArray();

                Points.AsParallel()
                    .Where(p => p.Tetrahedrons.All(t => t.IsIn[domaiNumber]) || p.Tetrahedrons.All(t => !t.IsIn[domaiNumber]))
                    .ForAll(p =>
                        {
                            foreach (int i in expressionIndexes)
                                p.Boundary[i] = false;
                        });

                Parallel.ForEach(Points, p =>
                {
                    foreach (int ineqNumber in p.Boundary.Cast<bool>().Select((b, i) => new { b = b, i = i }).Where(bi => bi.b).Select(bi => bi.i).ToArray())
                    {
                        if (!p.Points.Any(p1 => p1.Boundary[ineqNumber]))
                            p.Boundary[ineqNumber] = false;
                    }
                });
            }
        }

        private void ResolveIneq(int ineqNumber, int domaiNumber)
        {
            if (OnProgress != null)
                OnProgress(0.5 + 0.5 * (double)(ineqNumber + 1) / (ineqTreeBoxed.ExpressionList.Count));

            Func<Point, bool> isBoundaryPoint = (p => p.Tetrahedrons.Any(t => t.Boundary[ineqNumber]));

            Points.AsParallel().Where(p => isBoundaryPoint(p)).ForAll(p =>
            {
                Eval(p, ineqNumber);
                if (p.Boundary[ineqNumber])
                    p.U = 0;
            });

            ResolveEdges(EdgesForBoundary(ineqNumber).Where(e => (e.P1.U * e.P2.U < 0 && e.P1.HasAnyBoundary(e.P2))), ineqNumber, true);

            foreach (Point p in Points.AsParallel().Where(p => p.U == 0 && p.BoundaryCount > 0 && isBoundaryPoint(p)))
            {
                foreach (int bi in p.Boundary.Cast<bool>().Select((b, i) => new { b = b, i = i }).Where(bi => bi.b).Select(bi => bi.i))
                {
                    if (p.Points.Where(p1 => p1.Boundary[bi] && isBoundaryPoint(p1)).All(p1 => p1.U <= 0))
                    {
                        p.U = -1;
                        CenterPoint(p, 100, true);
                    }
                    else if (p.Points.Where(p1 => p1.Boundary[bi] && isBoundaryPoint(p1)).All(p1 => p1.U >= 0))
                    {
                        p.U = 1;
                        CenterPoint(p, 100, true);
                    }
                }
            }

            ResolveEdges(EdgesForBoundary(ineqNumber).Where(e => (e.P1.U * e.P2.U < 0)), ineqNumber, false);

            Edge[] innerEdges = Points.Where(p => p.U == 0 && isBoundaryPoint(p))
                    .SelectMany(
                        p1 => p1.Tetrahedrons.Where(t => t.Boundary[ineqNumber]).SelectMany(t => t.Points.Where(p2 => p2 > p1 && p2.U == 0)).Distinct(),
                        (p1, p2) => new Edge(p1, p2)
                        ).ToArray();

            foreach (Edge e in innerEdges)
            {
                if (!e.Valid)
                    continue;

                if (e.P1.Tetrahedrons.Where(t => t.Boundary[ineqNumber]).Intersect(e.P2.Tetrahedrons).All(t => t.Points.All(p => p.U >= 0)))
                {
                    Point newPoint = DivideEdge(e, -1, (e.P1 + e.P2) / 2);
                    newPoint.U = 1;
                    newPoint.Boundary[ineqNumber] = false;
                    CenterPoint(newPoint, 100, true);
                }
            }

            foreach (Edge e in innerEdges)
            {
                if (!e.Valid)
                    continue;

                if (e.P1.Tetrahedrons.Where(t => t.Boundary[ineqNumber]).Intersect(e.P2.Tetrahedrons).All(t => t.Points.All(p => p.U <= 0)))
                {
                    Point newPoint = DivideEdge(e, -1, (e.P1 + e.P2) / 2);
                    newPoint.U = -1;
                    newPoint.Boundary[ineqNumber] = false;
                    CenterPoint(newPoint, 100, true);
                }
            }

            foreach (Point p in Points.Where(p => isBoundaryPoint(p)))
            {
                if (p.U < 0)
                {
                    p.Boundary[ineqNumber] = false;
                }
                else if (p.U > 0)
                {
                    p.Boundary[ineqNumber] = false;
                }
                else
                {
                    if (p.Points.Where(p1 => isBoundaryPoint(p1)).All(p1 => p1.U <= 0))
                    {
                        p.U = -1;
                        p.Boundary[ineqNumber] = false;
                        CenterPoint(p, 100, true);
                    }
                    else if (p.Points.Where(p1 => isBoundaryPoint(p1)).All(p1 => p1.U >= 0))
                    {
                        p.U = 1;
                        p.Boundary[ineqNumber] = false;
                        CenterPoint(p, 100, true);
                    }
                    else
                    {
                        p.Boundary[ineqNumber] = true;
                    }
                }
            }

            Tetrahedrons.AsParallel().Where(t => t.Boundary[ineqNumber]).ForAll(t =>
            {
                t.IsIn[domaiNumber] = t.Points.Any(p => p.U < 0);// || t.Points.All(p => p.U == 0);
            });
        }

        private void ResolveEdges(ParallelQuery<Edge> qEdges, int ineqNumber, bool onSurface)
        {
            var edges = qEdges
                .Select(e =>
                            {
                                Point p = Bisection(e, ineqNumber);
                                Point nearPoint, farPoint;
                                bool movable;
                                double dist1, dist2, dist, length;

                                dist1 = p.Distance(e.P1);
                                dist2 = p.Distance(e.P2);
                                length = e.P1.Distance(e.P2);

                                if (dist1 < dist2)
                                {
                                    nearPoint = e.P1;
                                    farPoint = e.P2;
                                    movable = e.P1.HasAllBoundary(e.P2);
                                    dist = (length != 0 ? dist1 / length : 0);
                                }
                                else
                                {
                                    nearPoint = e.P2;
                                    farPoint = e.P1;
                                    movable = e.P2.HasAllBoundary(e.P1);
                                    dist = (length != 0 ? dist2 / length : 0);
                                }

                                if (dist >= maxWarpDist || (!nearPoint.Movable && dist >= 0.01d) || (nearPoint.Forced))
                                    movable = false;

                                if (nearPoint.BoundaryCount >= 3 && dist < 0.05) //less then 1/2^n in bisection is enough
                                    nearPoint.U = 0;

                                return new
                                {
                                    Edge = e,
                                    MidPoint = p,
                                    Dist = dist,
                                    NearPoint = nearPoint,
                                    FarPoint = farPoint,
                                    Movable = movable //&& nearPoint.Movable
                                };
                            }
                    )
                .ToArray();

            foreach (var e in edges.Where(e => e.Movable).OrderBy(e => e.Dist))
            {
                ResolveEdge(e.Edge, e.MidPoint, e.NearPoint, true, ineqNumber);
            }

            /*if (!onSurface)
            {
                foreach (var fanEdges in edges.Where(e => e.Edge.P1.U * e.Edge.P2.U < 0 && !e.Movable && e.NearPoint.BoundaryCount == 1 && e.Dist < maxWarpDist).GroupBy(e => e.NearPoint))
                {
                    var e = fanEdges.Where(e1 => e1.Edge.Valid).OrderBy(e1 => e1.Dist).FirstOrDefault();
                    
                    if (e == null)
                        continue;

                    ReplicateSurfacePoint(
                        e.NearPoint, 
                        e.FarPoint, 
                        fanEdges.Select(e1 => e1.Edge).ToArray(),
                        fanEdges.Select(e1 => e1.MidPoint).ToArray(),
                        e.NearPoint.BoundaryFirstIndex);                        
                }
            }*/

            foreach (var e in edges.Where(e => e.Edge.P1.U * e.Edge.P2.U < 0 && e.Edge.Valid))
            {
                ResolveEdge(e.Edge, e.MidPoint, e.NearPoint, false, ineqNumber);
            }
        }

        /*private Point ReplicateSurfacePoint(Point p, Point p1, Edge[] edges, Point[] midPoints,  int ineqNumber)
        {
            if (edges.Length < 2)
                return null;
            
            HashSet<Tetrahedron> fan = new HashSet<Tetrahedron>();

            foreach (var t in p.Tetrahedrons.Intersect(p1.Tetrahedrons))
                fan.Add(t);

            bool added = true;
            while(added)
            {
                added = false;
                foreach (var pp in fan.SelectMany(tt => tt.Points).Where(pp => pp != p && !pp.Boundary[ineqNumber]).ToArray())
                {
                    foreach (var t in p.Tetrahedrons.Intersect(pp.Tetrahedrons))
                    {
                        if(!fan.Contains(t))
                        {
                            fan.Add(t);
                            added = true;
                        }
                    }                    
                }
            }

            Point newPoint = new Point(p.X, p.Y, p.Z, p.Boundary.Length);
            newPoint.U = p.U;
            for (int i = 0; i < newPoint.Boundary.Length; i++)
            {
                newPoint.Boundary[i] = (i != ineqNumber) && p.Boundary[i];
            }

            AddPoint(newPoint);

            List<Tetrahedron> newTetrahedrons = new List<Tetrahedron>();

            foreach(var t in fan.Where(t => t.Points.Count(pp => pp != p && pp.Boundary[ineqNumber]) == 2))
            {
                Point[] fp = t.Points.Where(pp => pp != p && pp.Boundary[ineqNumber]).ToArray();

                newPoint.MoveTo(t.Points.Where(pp => !pp.Boundary[ineqNumber]).Single(), false);

                Tetrahedron tt = AddTetrahedron(newPoint, p, fp[0], fp[1], Math.Sign(t.Volume));
                tt.IsIn.Or(t.IsIn);
                tt.Boundary.Or(t.Boundary);
                tt.IsOnBoundary.Or(t.IsOnBoundary);               

                newTetrahedrons.Add(tt);
            }
            newPoint.MoveTo(p, false);

            foreach (var t in fan)
            {
                t.ReplacePoint(p, newPoint);
            }

            bool resolved = false;
            for (int i = 0; i < edges.Length; i++)
            {
                if (newPoint.CanMoveTo(midPoints[i]))
                {
                    ResolveEdge(edges[i], midPoints[i], newPoint, true, ineqNumber);
                    resolved = true;
                }
            }

            if (!resolved)
            {
                foreach (var t in fan)
                {
                    t.ReplacePoint(newPoint, p);
                }
                foreach (Tetrahedron t in newTetrahedrons)
                {
                    DeleteTetrahedron(t);
                }
            }

            return newPoint;
        }*/

        private void ResolveEdge(Edge e, Point midPoint, Point nearPoint, bool movable, int ineqNumber)
        {
            if (e.P1.U * e.P2.U >= 0)
                return;

            if (nearPoint.Tetrahedrons.Any(t => t.Points.Where(p => p.U == 0).Count() >= 3))
            {
                movable = false;
            }

            if (movable)
            {
                movable = nearPoint.MoveTo(midPoint, true);
                if (movable)
                {
                    nearPoint.U = 0;

                    if (nearPoint.BoundaryCount == 1)
                        ProjectToEdge(nearPoint, ineqNumber, nearPoint.BoundaryFirstIndex, true);
                    else if (nearPoint.BoundaryCount == 2)
                    {
                        ProjectToCorner(nearPoint, ineqNumber, nearPoint.BoundaryFirstIndex, nearPoint.BoundarySecondIndex, true);
                    }
                }
            }
            if (!movable)
            {
                if (Math.Abs(e.P1.U) < 1e-10)
                {
                    e.P1.U = 0;
                    return;
                }
                if (Math.Abs(e.P2.U) < 1e-10)
                {
                    e.P2.U = 0;
                    return;
                }
                DivideEdge(e, ineqNumber, midPoint);
            }
        }

        public Point DivideEdge(Edge e, int ineqNumber, Point newPoint)
        {
            if (newPoint == null)
                newPoint = Bisection(e, ineqNumber);

            for (int i = 0; i < newPoint.Boundary.Length; i++)
            {
                newPoint.Boundary[i] = e.P1.Boundary[i] && e.P2.Boundary[i];
            }

            newPoint.Movable = e.P1.Movable && e.P2.Movable;
            newPoint.Forced = e.P1.Forced && e.P2.Forced;

            newPoint.U = 0;
            AddPoint(newPoint);

            foreach (Tetrahedron t in e.P1.Tetrahedrons.Intersect(e.P2.Tetrahedrons).ToArray())
            {
                Point[] p = t.Points.Where(pp => pp != e.P1 && pp != e.P2).ToArray();

                int volumeSign = Math.Sign(t.Volume);

                Tetrahedron tt = AddTetrahedron(p[0], p[1], e.P1, newPoint, volumeSign);
                tt.IsIn.Or(t.IsIn);
                tt.Boundary.Or(t.Boundary);
                tt.IsOnBoundary.Or(t.IsOnBoundary);


                tt = AddTetrahedron(p[0], p[1], e.P2, newPoint, volumeSign);
                tt.IsIn.Or(t.IsIn);
                tt.Boundary.Or(t.Boundary);
                tt.IsOnBoundary.Or(t.IsOnBoundary);


                DeleteTetrahedron(t);
            }

            if (newPoint.BoundaryCount == 1 && ineqNumber >= 0)
                ProjectToEdge(newPoint, ineqNumber, newPoint.BoundaryFirstIndex, true);
            else if (newPoint.BoundaryCount == 2 && ineqNumber >= 0)
            {
                ProjectToCorner(newPoint, ineqNumber, newPoint.BoundaryFirstIndex, newPoint.BoundarySecondIndex, true);
            }

            return newPoint;
        }

        public bool CollapseEdge(Edge e, Point toPoint, bool jiggle = false)
        {
            Point fromPoint = null;

            if (toPoint == null)
            {
                if (e.P2.HasAllBoundary(e.P1))
                {
                    toPoint = e.P1;
                    fromPoint = e.P2;
                }
                else if (e.P1.HasAllBoundary(e.P2))
                {
                    toPoint = e.P2;
                    fromPoint = e.P1;
                }
                else
                {
                    return false;
                }
            }
            else
            {
                fromPoint = (toPoint == e.P1 ? e.P2 : e.P1);
                if (!fromPoint.HasAllBoundary(toPoint))
                    return false;
            }

            if (!fromPoint.CanMoveTo(toPoint, fromPoint.Tetrahedrons.Intersect(toPoint.Tetrahedrons), 0.15d))
            {
                return false;
            }

            foreach (Tetrahedron t in fromPoint.Tetrahedrons.Intersect(toPoint.Tetrahedrons).ToArray())
            {
                DeleteTetrahedron(t);
            }

            foreach (Tetrahedron t in fromPoint.Tetrahedrons.ToArray())
            {
                t.ReplacePoint(fromPoint, toPoint);
            }

            DeletePoint(fromPoint);

            if (jiggle)
            {
                Jiggle(2, toPoint.Points.SelectMany(p => p.Points));
            }
            return true;
        }

        protected Point Bisection(Edge e, int ineqNumber)
        {
            Point a = new Point(e.P1.X, e.P1.Y, e.P1.Z, e.P1.Boundary.Length);
            Point b = new Point(e.P2.X, e.P2.Y, e.P2.Z, e.P1.Boundary.Length);
            Point c = null;
            double ua = e.P1.U;
            double ub = e.P2.U;
            double uc = 0;

            for (int i = 0; i < 8; i++)
            {
                c = (a + b) / 2;
                uc = Eval(c.X, c.Y, c.Z, ineqNumber);

                if (uc == 0)
                    return c;

                if (ua * uc < 0)
                {
                    ub = uc;
                    b = c;
                }
                else
                {
                    ua = uc;
                    a = c;
                }
            }

            return c;
        }

        public bool ProjectToSurface(Point P, double precision, int ineqNumber, bool safe)
        {
            if (ProjectToSurfaceSpec.Count > 0 && ProjectToSurfaceSpec.ContainsKey(ineqTreeBoxed.ExpressionList[ineqNumber]))
            {
                Point P1 = new Point(P.X, P.Y, P.Z);

                if (ProjectToSurfaceSpec[ineqTreeBoxed.ExpressionList[ineqNumber]](P1))
                {
                    return P.MoveTo(P1, safe);
                }
            }

            double n1, n2, n3, n;
            double dx, dy, dz;
            double w, wx, wy, wz;
            dx = D / 10000000d;
            dy = D / 10000000d;
            dz = D / 10000000d;

            w = Eval(P.X, P.Y, P.Z, ineqNumber);

            wx = Eval(P.X + dx, P.Y, P.Z, ineqNumber);
            wy = Eval(P.X, P.Y + dy, P.Z, ineqNumber);
            wz = Eval(P.X, P.Y, P.Z + dz, ineqNumber);
            n1 = (wx - w) / dx;
            n2 = (wy - w) / dy;
            n3 = (wz - w) / dz;

            n = Math.Sqrt(n1 * n1 + n2 * n2 + n3 * n3);
            if (n != 0)
            {
                n1 = n1 / n; n2 = n2 / n; n3 = n3 / n;
            }
            else
            {
                return false;
            }
            Point A, B;
            double w1, tet;

            A = new Point(P.X, P.Y, P.Z);
            tet = -w / n;
            B = new Point(P.X + tet * n1, P.Y + tet * n2, P.Z + tet * n3);

            double stopdist = D / precision;
            double errordist = 2 * D / 3.0d;
            int p = 0;
            int maxItCount = 100;

            if (tet > errordist)
            {
                p = maxItCount;
            }

            while (Math.Abs(tet) > stopdist && p < maxItCount)
            {
                p++;
                w1 = Eval(B.X, B.Y, B.Z, ineqNumber);
                if (w1 == w || tet > errordist)
                {
                    p = maxItCount;
                    break;
                }
                tet = -w1 * tet / (w1 - w);
                A.MoveTo(B, false);
                w = w1;
                B.MoveTo(B.X + tet * n1, B.Y + tet * n2, B.Z + tet * n3, false);
            }
            if (p < maxItCount && Math.Sqrt((P.X - B.X) * (P.X - B.X) + (P.Y - B.Y) * (P.Y - B.Y) + (P.Z - B.Z) * (P.Z - B.Z)) < errordist)
            {
                return P.MoveTo(B, safe);
            }
            else
                return false;
        }

        public bool ProjectToEdge(Point P, int ineqNumber1, int ineqNumber2, bool safe)
        {
            double nx1, ny1, nz1;//, n1;
            double nx2, ny2, nz2;//, n2;
            double dx, dy, dz;
            double w, wx, wy, wz;
            dx = D / 100000d;
            dy = D / 100000d;
            dz = D / 100000d;

            w = Eval(P.X, P.Y, P.Z, ineqNumber1);
            wx = Eval(P.X + dx, P.Y, P.Z, ineqNumber1);
            wy = Eval(P.X, P.Y + dy, P.Z, ineqNumber1);
            wz = Eval(P.X, P.Y, P.Z + dz, ineqNumber1);
            nx1 = (wx - w) / dx;
            ny1 = (wy - w) / dy;
            nz1 = (wz - w) / dz;
            //n1 = Math.Sqrt(nx1 * nx1 + ny1 * ny1 + nz1 * nz1);

            w = Eval(P.X, P.Y, P.Z, ineqNumber2);
            wx = Eval(P.X + dx, P.Y, P.Z, ineqNumber2);
            wy = Eval(P.X, P.Y + dy, P.Z, ineqNumber2);
            wz = Eval(P.X, P.Y, P.Z + dz, ineqNumber2);
            nx2 = (wx - w) / dx;
            ny2 = (wy - w) / dy;
            nz2 = (wz - w) / dz;
            //n2 = Math.Sqrt(nx1 * nx1 + ny1 * ny1 + nz1 * nz1);

            double hh1 = 0, hh2 = 0;
            if (Numeric.Newton2d(ref hh1, ref hh2,
                (h1, h2) => Eval(P.X + h1 * nx1 + h2 * nx2, P.Y + h1 * ny1 + h2 * ny2, P.Z + h1 * nz1 + h2 * nz2, ineqNumber1),
                (h1, h2) => Eval(P.X + h1 * nx1 + h2 * nx2, P.Y + h1 * ny1 + h2 * ny2, P.Z + h1 * nz1 + h2 * nz2, ineqNumber2), D))
            {
                //******** safe max D distance move test added ***************
                if (Math.Sqrt((hh1 * nx1 + hh2 * nx2) * (hh1 * nx1 + hh2 * nx2) +
                (hh1 * ny1 + hh2 * ny2) * (hh1 * ny1 + hh2 * ny2) +
                (hh1 * nz1 + hh2 * nz2) * (hh1 * nz1 + hh2 * nz2)) < D)
                {

                    return P.MoveTo(
                            P.X + hh1 * nx1 + hh2 * nx2,
                            P.Y + hh1 * ny1 + hh2 * ny2,
                            P.Z + hh1 * nz1 + hh2 * nz2,
                            safe
                        );
                }
                else
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
        }

        protected bool ProjectToCorner(Point P, int ineqNumber1, int ineqNumber2, int ineqNumber3, bool safe)
        {
            double x0 = P.X, y0 = P.Y, z0 = P.Z;
            if (Numeric.Newton3d(ref x0, ref y0, ref z0, (x, y, z) => Eval(x, y, z, ineqNumber1), (x, y, z) => Eval(x, y, z, ineqNumber2), (x, y, z) => Eval(x, y, z, ineqNumber3), D))
            {
                return P.MoveTo(x0, y0, z0, safe);
            }
            else
                return false;
        }

        private void CreateBackgroundMesh()
        {
            double x0, y0, z0, x1, y1, z1;

            x0 = X0 - D;
            x1 = X1 + D;
            y0 = Y0 - D;
            y1 = Y1 + D;
            z0 = Z0 - D;
            z1 = Z1 + D;

            int xCount = (int)Math.Ceiling(((x1 - x0) / D));
            int yCount = (int)Math.Ceiling(((y1 - y0) / D));
            int zCount = (int)Math.Ceiling(((z1 - z0) / D));

            Point[,,] points = new Point[xCount + 1, yCount + 1, zCount + 1];
            Point[,,] shiftPoints = new Point[xCount + 1, yCount + 1, zCount + 1];

            for (int k = 0; k <= zCount; k++)
            {
                for (int j = 0; j <= yCount; j++)
                {
                    for (int i = 0; i <= xCount; i++)
                    {
                        points[i, j, k] = CreatePoint(x0 + i * D, y0 + j * D, z0 + k * D);
                        shiftPoints[i, j, k] = CreatePoint(x0 + i * D + D / 2.0d, y0 + j * D + D / 2.0d, z0 + k * D + D / 2.0d);
                    }
                }
            }

            this.points = new HashSet<Point>(points.Cast<Point>().Concat(shiftPoints.Cast<Point>()));

            Tetrahedron[] tetrahedrons = new Tetrahedron[xCount * yCount * zCount * 12];
            int p = 0;
            for (int k = 0; k < zCount; k++)
            {
                for (int j = 0; j < yCount; j++)
                {
                    for (int i = 0; i < xCount; i++)
                    {
                        tetrahedrons[p++] = CreateTetrahedron(shiftPoints[i, j, k], shiftPoints[i + 1, j, k], points[i + 1, j, k], points[i + 1, j + 1, k]);
                        tetrahedrons[p++] = CreateTetrahedron(shiftPoints[i, j, k], shiftPoints[i + 1, j, k], points[i + 1, j, k + 1], points[i + 1, j + 1, k + 1]);
                        tetrahedrons[p++] = CreateTetrahedron(shiftPoints[i, j, k], shiftPoints[i + 1, j, k], points[i + 1, j, k], points[i + 1, j, k + 1]);
                        tetrahedrons[p++] = CreateTetrahedron(shiftPoints[i, j, k], shiftPoints[i + 1, j, k], points[i + 1, j + 1, k], points[i + 1, j + 1, k + 1]);

                        tetrahedrons[p++] = CreateTetrahedron(shiftPoints[i, j, k], shiftPoints[i, j + 1, k], points[i, j + 1, k], points[i + 1, j + 1, k]);
                        tetrahedrons[p++] = CreateTetrahedron(shiftPoints[i, j, k], shiftPoints[i, j + 1, k], points[i, j + 1, k + 1], points[i + 1, j + 1, k + 1]);
                        tetrahedrons[p++] = CreateTetrahedron(shiftPoints[i, j, k], shiftPoints[i, j + 1, k], points[i, j + 1, k], points[i, j + 1, k + 1]);
                        tetrahedrons[p++] = CreateTetrahedron(shiftPoints[i, j, k], shiftPoints[i, j + 1, k], points[i + 1, j + 1, k], points[i + 1, j + 1, k + 1]);

                        tetrahedrons[p++] = CreateTetrahedron(shiftPoints[i, j, k], shiftPoints[i, j, k + 1], points[i, j, k + 1], points[i, j + 1, k + 1]);
                        tetrahedrons[p++] = CreateTetrahedron(shiftPoints[i, j, k], shiftPoints[i, j, k + 1], points[i + 1, j, k + 1], points[i + 1, j + 1, k + 1]);
                        tetrahedrons[p++] = CreateTetrahedron(shiftPoints[i, j, k], shiftPoints[i, j, k + 1], points[i, j, k + 1], points[i + 1, j, k + 1]);
                        tetrahedrons[p++] = CreateTetrahedron(shiftPoints[i, j, k], shiftPoints[i, j, k + 1], points[i, j + 1, k + 1], points[i + 1, j + 1, k + 1]);
                    }
                }
            }
            this.tetrahedrons = new HashSet<Tetrahedron>(tetrahedrons);
        }

        public IEnumerable<Point> RefineTetrahedralMeshRedGreen(IEnumerable<Tetrahedron> refineList)
        {
            var newPoints = new Dictionary<Edge, Point>();
            var refine = refineList.ToHashSet();

            foreach (var t in refineList.ToArray())
            {
                RefineRed(t, newPoints, refine);
            }

            bool red = true;
            while (red)
            {
                red = false;
                foreach (var t in refine.ToArray())
                {
                    var dividedEdges = t.Edges().Where(e => newPoints.ContainsKey(e)).ToArray();

                    if (!(
                            dividedEdges.Length == 1 ||
                            (dividedEdges.Length == 2 && !dividedEdges[0].Intersect(dividedEdges[1]).Any()) ||
                            (dividedEdges.Length == 3 && dividedEdges.SelectMany(e => e).Distinct().Count() == 3)
                        )
                    )
                    {
                        red = true;
                        RefineRed(t, newPoints, refine);
                    }
                }
            }

            foreach (var t in refine/*.ToArray()*/)
            {
                var dividedEdges = t.Edges().Where(e => newPoints.ContainsKey(e)).ToArray();

                int volumeSign = Math.Sign(t.Volume);
                var newTetras = new List<Tetrahedron>();

                if (dividedEdges.Length == 1)
                {
                    var e = dividedEdges[0];
                    var np = newPoints[e];
                    var v = t.Points.Where(p => p != e.P1 && p != e.P2).ToArray();

                    newTetras.Add(AddTetrahedron(np, e.P1, v[0], v[1], volumeSign));
                    newTetras.Add(AddTetrahedron(np, e.P2, v[0], v[1], volumeSign));
                }

                if (dividedEdges.Length == 2 && !dividedEdges[0].Intersect(dividedEdges[1]).Any())
                {
                    var e1 = dividedEdges[0];
                    var e2 = dividedEdges[1];
                    var p1 = newPoints[e1];
                    var p2 = newPoints[e2];

                    newTetras.Add(AddTetrahedron(e1.P1, p1, e2.P1, p2, volumeSign));
                    newTetras.Add(AddTetrahedron(e1.P2, p1, e2.P1, p2, volumeSign));
                    newTetras.Add(AddTetrahedron(e1.P1, p1, e2.P2, p2, volumeSign));
                    newTetras.Add(AddTetrahedron(e1.P2, p1, e2.P2, p2, volumeSign));
                }

                if (dividedEdges.Length == 3 && dividedEdges.SelectMany(e => e).Distinct().Count() == 3)
                {
                    var p0 = t.Points.Where(p => !dividedEdges.Any(e => e.Contains(p))).Single();
                    var mp = dividedEdges.Select(e => newPoints[e]).ToArray();

                    newTetras.Add(AddTetrahedron(p0, mp[0], mp[1], mp[2], volumeSign));
                    newTetras.Add(AddTetrahedron(p0, mp[0], mp[1], dividedEdges[0].Intersect(dividedEdges[1]).Single(), volumeSign));
                    newTetras.Add(AddTetrahedron(p0, mp[0], mp[2], dividedEdges[0].Intersect(dividedEdges[2]).Single(), volumeSign));
                    newTetras.Add(AddTetrahedron(p0, mp[1], mp[2], dividedEdges[1].Intersect(dividedEdges[2]).Single(), volumeSign));
                }

                foreach (var tt in newTetras)
                {
                    tt.IsIn.Or(t.IsIn);
                    tt.Boundary.Or(t.Boundary);
                    tt.IsOnBoundary.Or(t.IsOnBoundary);
                }

                DeleteTetrahedron(t);
                /*if (refine.Contains(t))
                {
                    refine.Remove(t);
                }*/
            }

            return newPoints.Values;
        }

        public void RefineRed(Tetrahedron t, Dictionary<Edge, Point> newPoints, HashSet<Tetrahedron> refine)
        {
            if (!Tetrahedrons.Contains(t))
            {
                return;
            }

            var midpoints = new Dictionary<(int, int), Point>();

            for (int i = 0; i < 4; i++)
            {
                for (int j = i + 1; j < 4; j++)
                {
                    var e = new Edge(t.Points[i], t.Points[j]);

                    if (newPoints.ContainsKey(e))
                    {
                        midpoints[(i, j)] = newPoints[e];
                    }
                    else
                    {
                        var np = (e.P1 + e.P2) / 2;
                        for (int b = 0; b < np.Boundary.Length; b++)
                        {
                            np.Boundary[b] = e.P1.Boundary[b] && e.P2.Boundary[b];
                        }
                        np.Movable = e.P1.Movable && e.P2.Movable;
                        np.Forced = e.P1.Forced && e.P2.Forced;
                        np.U = 0;
                        AddPoint(np);
                        midpoints[(i, j)] = np;
                        newPoints[e] = np;

                        foreach (Tetrahedron tt in e.P1.Tetrahedrons.Intersect(e.P2.Tetrahedrons).ToArray())
                        {
                            if (!refine.Contains(tt))
                            {
                                refine.Add(tt);
                            }
                        }
                    }
                }
            }

            int volumeSign = Math.Sign(t.Volume);
            var newTetras = new List<Tetrahedron>();
            var v = t.Points;
            var m = midpoints;

            newTetras.Add(AddTetrahedron(v[0], m[(0, 1)], m[(0, 2)], m[(0, 3)], volumeSign));
            newTetras.Add(AddTetrahedron(v[1], m[(0, 1)], m[(1, 2)], m[(1, 3)], volumeSign));
            newTetras.Add(AddTetrahedron(v[2], m[(0, 2)], m[(1, 2)], m[(2, 3)], volumeSign));
            newTetras.Add(AddTetrahedron(v[3], m[(0, 3)], m[(1, 3)], m[(2, 3)], volumeSign));

            newTetras.Add(AddTetrahedron(m[(0, 1)], m[(0, 2)], m[(0, 3)], m[(1, 2)], volumeSign));
            newTetras.Add(AddTetrahedron(m[(1, 2)], m[(1, 3)], m[(2, 3)], m[(0, 3)], volumeSign));
            newTetras.Add(AddTetrahedron(m[(0, 2)], m[(1, 2)], m[(2, 3)], m[(0, 3)], volumeSign));
            newTetras.Add(AddTetrahedron(m[(0, 1)], m[(1, 2)], m[(1, 3)], m[(0, 3)], volumeSign));

            foreach (var tt in newTetras)
            {
                tt.IsIn.Or(t.IsIn);
                tt.Boundary.Or(t.Boundary);
                tt.IsOnBoundary.Or(t.IsOnBoundary);
            }

            DeleteTetrahedron(t);
            if (refine.Contains(t))
            {
                refine.Remove(t);
            }
        }

        //public void RefineTetrahedralMeshNearPoint(Point c, double r)
        //{
        //    var tetras = Tetrahedrons.Where(t => /*t.IsOnBoundaryDomain(0) && */t.Points.Any(p => p.Distance(c) <= r)).ToArray();

        //    foreach (var t in tetras)
        //    {
        //        if (!Tetrahedrons.Contains(t))
        //            continue;

        //        var ee = t.Edges().OrderByDescending(e => e.SqrLength).First();

        //        var p = DivideEdge(ee, -1, (ee.P1 + ee.P2) / 2);
        //    }
        //}

        public void RefineTetrahedralMeshNearPointRedGreen(IEnumerable<Point> c, double r, int count = 1)
        {
            for (int i = 0; i < count; i++)
            {
                var tetras = Tetrahedrons.Where(t => c.Any(p => t.Points.Average().Distance(p) <= r)).ToArray();
                RefineTetrahedralMeshRedGreen(tetras);
            }
        }

        //public void RefineTetrahedralMeshByTetrahedrons(int ineqNumber, int count = 5, double tolerance = 0.1d, bool redgreen = false)
        //{
        //    Points.AsParallel().ForAll(p => p.U = 0);

        //    Parallel.ForEach(Points.Where(p => p.Tetrahedrons.Any(t => t.Boundary[ineqNumber])), p =>
        //    {
        //        Eval(p, ineqNumber);
        //    });

        //    var good = new HashSet<Tetrahedron>();


        //    for (int i = 0; i < count; i++)
        //    {
        //        var tetras = Tetrahedrons.AsParallel().Where(t =>
        //        {
        //            if (
        //                !t.Boundary[ineqNumber] || good.Contains(t) 
        //                || 
        //                (
        //                    //redgreen && i < 1 &&
        //                    (
        //                        t.Points.All(p =>  Math.Sign(p.U) == 1) || t.Points.All(p => Math.Sign(p.U) == -1)
        //                    )
        //                ))
        //            {
        //                return false;
        //            }

        //            var midP = t.Points.Average();
        //            Eval(midP, ineqNumber);

        //            var res = Math.Abs(midP.U - t.Points.Sum(p => p.U) / 4.0d) > tolerance * D; /*t.Edges().Max(e => Math.Sqrt(e.SqrLength))*/ /*(t.Points.Max(p => p.U) - t.Points.Min(p => p.U))*/

        //            if (!res)
        //            {
        //                lock (good)
        //                {
        //                    good.Add(t);
        //                }
        //            }

        //            return res;

        //        })
        //        .Select(t => 
        //        {
        //            var ee = t.Edges().OrderByDescending(e => e.SqrLength).First();
        //            return new 
        //            { 
        //                t,
        //                ee
        //            }; 
        //        })
        //        .OrderByDescending(t => t.ee.SqrLength)
        //        .ToArray();

        //        if (tetras.Length == 0)
        //        {
        //            break;
        //        }

        //        if (redgreen && i < 1)
        //        {
        //            var newPoints = RefineTetrahedralMeshRedGreen(tetras.Select(t => t.t));
        //            Parallel.ForEach(newPoints, p =>
        //            {
        //                Eval(p, ineqNumber);
        //                p.Movable = false;
        //                foreach (var pp in p.Points)
        //                {
        //                    pp.Movable = false;
        //                }
        //            });
        //        }
        //        else
        //        {
        //            foreach (var t in tetras)
        //            {
        //                if (!Tetrahedrons.Contains(t.t))
        //                    continue;

        //                var ee = t.ee;

        //                var p = DivideEdge(ee, -1, (ee.P1 + ee.P2) / 2);
        //                Eval(p, ineqNumber);
        //                foreach (var tt in p.Tetrahedrons) 
        //                {
        //                    tt.Boundary[ineqNumber] = true;
        //                }
        //                foreach(var pp in p.Points)
        //                {
        //                    if (pp.U == 0)
        //                    {
        //                        Eval(pp, ineqNumber);
        //                    }
        //                }

        //                p.Movable = false;
        //                foreach (var pp in p.Points)
        //                {
        //                    pp.Movable = false;
        //                }
        //            }
        //        }
        //    }
        //}


        public void RefineTetrahedralMesh(IneqTree ineq, int count = 5, double tolerance = 0.1d, bool redgreen = false, double relativeTolerance = 0)
        {
            RefineTetrahedralMesh(ineq.Root, count, tolerance, redgreen, relativeTolerance);
        }

        private void RefineTetrahedralMesh(IneqTree.IneqNode ineq, int count, double tolerance, bool redgreen, double relativeTolerance)
        {
            if (ineq.NodeType == IneqTree.NodeType.NodeExpression)
            {
                int ineqNumber = ineqTreeBoxed.ExpressionList.IndexOf(ineq.Expression);
                if (ineqNumber >= 0)
                {
                    RefineTetrahedralMesh(ineqNumber, count, tolerance, redgreen, relativeTolerance);
                }
            }
            else
            {
                RefineTetrahedralMesh(ineq.Left, count, tolerance, redgreen, relativeTolerance);
                RefineTetrahedralMesh(ineq.Right, count, tolerance, redgreen, relativeTolerance);
            }
        }

        public void RefineTetrahedralMesh(int ineqNumber, int count = 5, double tolerance = 0.1d, bool redgreen = false, double relativeTolerance = 0)
        {
            Parallel.ForEach(Points, p =>
            {
                Eval(p, ineqNumber);
            });

            var good = new HashSet<Tetrahedron>();

            for (int i = 0; i < count; i++)
            {
                var tetras = Tetrahedrons.AsParallel().Where(t =>
                {
                    if (!t.IntersectBoundary() && !t.CanIntersectBoundary())
                    {
                        return false;
                    }

                    var midP = t.Points.Average();
                    Eval(midP, ineqNumber);

                    var err = Math.Abs(midP.U - t.Points.Sum(p => p.U) / 4.0d);
                    var res = (tolerance > 0 && err > tolerance * D) || (relativeTolerance > 0 && err > relativeTolerance * Math.Sqrt(t.Edges().Max(e => e.SqrLength)));

                    if (!res)
                    {
                        lock (good)
                        {
                            good.Add(t);
                        }
                    }

                    return res;
                })
                .Select(t =>
                {
                    var ee = t.Edges().OrderByDescending(e => e.SqrLength).First();
                    return new
                    {
                        t,
                        ee
                    };
                })
                .OrderByDescending(t => t.ee.SqrLength)
                .ToArray();

                if (tetras.Length == 0)
                {
                    break;
                }

                if (redgreen && i < 1)
                {
                    var newPoints = RefineTetrahedralMeshRedGreen(tetras.Select(t => t.t));
                    Parallel.ForEach(newPoints, p =>
                    {
                        Eval(p, ineqNumber);
                        p.Movable = false;
                        foreach (var pp in p.Points)
                        {
                            pp.Movable = false;
                        }
                    });
                }
                else
                {
                    foreach (var t in tetras)
                    {
                        if (!Tetrahedrons.Contains(t.t))
                            continue;

                        var ee = t.ee;

                        var newPoint = (ee.P1 + ee.P2) / 2;
                        double minQuality = 1.0d;
                        foreach (Tetrahedron tt in ee.P1.Tetrahedrons.Intersect(ee.P2.Tetrahedrons).ToArray())
                        {
                            Point[] pp = tt.Points.Where(ppp => ppp != ee.P1 && ppp != ee.P2).ToArray();

                            int volumeSign = Math.Sign(tt.Volume);

                            Tetrahedron ttt = new Tetrahedron(pp[0], pp[1], ee.P1, newPoint, 0, 0, true, volumeSign);
                            minQuality = Math.Min(minQuality, ttt.Quality);
                            ttt = new Tetrahedron(pp[0], pp[1], ee.P2, newPoint, 0, 0, true, volumeSign);
                            minQuality = Math.Min(minQuality, ttt.Quality);
                        }

                        if (minQuality < 0.01d)
                        {
                            continue;
                        }

                        var p = DivideEdge(ee, -1, (ee.P1 + ee.P2) / 2);
                        Eval(p, ineqNumber);

                        p.Movable = false;
                        foreach (var pp in p.Points)
                        {
                            pp.Movable = false;
                        }
                    }
                }
            }
        }

        public void FindCrossEdges(int ineqNumber, int division = 2)
        {
            Parallel.ForEach(Points.Where(p => p.Tetrahedrons.Any(t => t.Boundary[ineqNumber])), p =>
            {
                Eval(p, ineqNumber);
            });

            int find = 1;

            //while (find > 0)
            {
                find = 0;

                var edges = EdgesForBoundary(ineqNumber).Where(e =>
                {
                    var res = e.P1.U * e.P2.U > 0;
                    return res;
                })
                .OrderByDescending(e => e.SqrLength)
                .ToArray();

                foreach (var e in edges)
                {
                    var l = Math.Sqrt(e.SqrLength);
                    var d1 = Math.Abs(e.P1.U);
                    var d2 = Math.Abs(e.P2.U);

                    if (d1 + d2 > l)
                    {
                        continue;
                    }

                    var p1 = e.P1 + (d1 / l) * (e.P2 - e.P1);
                    var p2 = e.P2 + (d2 / l) * (e.P1 - e.P2);

                    for (int i = 1; i < division; i++)
                    {
                        var p = p1 + ((double)i / (double)division) * (p2 - p1);

                        Eval(p, ineqNumber);

                        if (e.P1.U * p.U < 0)
                        {
                            find++;
                            var u = p.U;
                            p = DivideEdge(e, -1, p);
                            p.U = u;
                            p.Tetrahedrons.AsParallel().ForAll(tt => tt.Boundary[ineqNumber] = true);

                            p.Movable = false;
                            foreach (var pp in p.Points)
                            {
                                pp.Movable = false;
                            }

                            break;
                        }

                    }
                }
            }
        }

        //public void RefineTetrahedralMeshByEdges(int ineqNumber, int count = 2, double tolerance = 0.1d)
        //{
        //    Parallel.ForEach(Points.Where(p => p.Tetrahedrons.Any(t => t.Boundary[ineqNumber])), p =>
        //    {
        //        Eval(p, ineqNumber);
        //    });

        //    var good = new HashSet<Edge>();

        //    for (int i = 0; i < count; i++)
        //    {
        //        var edges = EdgesForBoundary(ineqNumber).Where(e =>
        //        {
        //            if (count > 1 && good.Contains(e))
        //            {
        //                return false;
        //            }

        //            var midP = (e.P1 + e.P2) / 2;
        //            Eval(midP, ineqNumber);

        //            var res = Math.Abs(midP.U - (e.P1.U + e.P2.U) / 2.0d) > tolerance * D; // Math.Sqrt(e.SqrLength);

        //            if (count > 1 && !res)
        //            {
        //                lock (good)
        //                {
        //                    good.Add(e);
        //                }
        //            }

        //            return res;

        //        })
        //        .OrderByDescending(e => e.SqrLength)
        //        .ToArray();

        //        foreach (var e in edges)
        //        {
        //            var p = DivideEdge(e, -1, (e.P1 + e.P2) / 2);
        //            Eval(p, ineqNumber);
        //            p.Tetrahedrons.AsParallel().ForAll(tt => tt.Boundary[ineqNumber] = true);
        //        }

        //    }
        //}

        //public void RefineTetrahedralMeshByEdgesMax(int ineqNumber, int count = 5, double tolerance = 0.1d)
        //{
        //    Parallel.ForEach(Points.Where(p => p.Tetrahedrons.Any(t => t.Boundary[ineqNumber])), p =>
        //    {
        //        Eval(p, ineqNumber);
        //    });

        //    var good = new HashSet<Edge>();

        //    for (int i = 0; i < count; i++)
        //    {
        //        var edges = EdgesForBoundary(ineqNumber).Where(e =>
        //        {
        //            if (count > 1 && good.Contains(e))
        //            {
        //                return false;
        //            }

        //            var midP = (e.P1 + e.P2) / 2;
        //            Eval(midP, ineqNumber);

        //            var res = Math.Abs(midP.U - (e.P1.U + e.P2.U) / 2.0d) > tolerance * D; /*Math.Sqrt(e.SqrLength)*/

        //            if (count > 1 && !res)
        //            {
        //                lock (good)
        //                {
        //                    good.Add(e);
        //                }
        //            }

        //            return res;

        //        })
        //        .Select(e =>
        //        {
        //            double tmax = Numeric.FindMaximum(t =>
        //            {
        //                Point pt = e.P1 + t * (e.P2 - e.P1);
        //                Eval(pt, ineqNumber);
        //                return Math.Abs(pt.U - (e.P1.U + t * (e.P2.U - e.P1.U)));
        //            });

        //            return new
        //            {
        //                e,
        //                tmax
        //            };
        //        })
        //        .OrderBy(e => e.tmax < 0.1d || e.tmax > 0.9d ? 0 : 1)
        //        .ThenByDescending(e => e.e.SqrLength)
        //        .ToArray();

        //        foreach (var ee in edges)
        //        {
        //            var e = ee.e;
        //            var tmax = ee.tmax;

        //            Point pmax = e.P1 + tmax * (e.P2 - e.P1);

        //            bool moved = false;

        //            if (tmax < 0.1d)
        //            {
        //                moved = e.P1.MoveTo(pmax, true);
        //                //e.P1.Movable = false;
        //            }
        //            else if (tmax > 0.9d)
        //            {
        //                moved = e.P2.MoveTo(pmax, true);
        //                //e.P2.Movable = false;
        //            }

        //            if(!moved)
        //            {
        //                var p = DivideEdge(e, -1, pmax);
        //                Eval(p, ineqNumber);
        //                //p.Movable = false;
        //                p.Tetrahedrons.AsParallel().ForAll(tt => tt.Boundary[ineqNumber] = true);
        //            }
        //        }
        //    }
        //}

        public void Jiggle(int count, bool edges = true)
        {
            if (Point.EnableUnsafeMove)
            {
                for (int i = 0; i < count; i++)
                    Points.AsParallel().ForAll(point =>
                    {
                        CenterPoint(point, i == count - 1 ? 1000 : 100, true, edges);
                    });
            }
            else
            {
                for (int i = 0; i < count; i++)
                    foreach (var point in Points)
                    {
                        CenterPoint(point, i == count - 1 ? 1000 : 100, true, edges);
                    }
            }
        }

        public void Jiggle(int count, IEnumerable<Point> points, bool edges = true)
        {
            if (Point.EnableUnsafeMove)
            {
                for (int i = 0; i < count; i++)
                {
                    points.AsParallel().ForAll(point =>
                    {
                        CenterPoint(point, i == count - 1 ? 1000 : 100, true, edges);
                    });
                }
            }
            else
            {
                for (int i = 0; i < count; i++)
                {
                    foreach (var point in points)
                    {
                        CenterPoint(point, i == count - 1 ? 1000 : 100, true, edges);
                    }
                }
            }
        }

        public void JiggleBackgroundMash(int count)
        {
            var boundaryTriangles = Tetrahedrons.SelectMany(t => t.Triangles().Where(tr => tr.Boundary))
                .ToArray();

            var boundaryPoints = boundaryTriangles.SelectMany(t => t).Distinct().ToHashSet();

            for (int i = 0; i < count; i++)
                Points.AsParallel().Where(p => !boundaryPoints.Contains(p) && !p.Forced).ForAll(point =>
                //foreach(var point in Points.Where(p => !boundaryPoints.Contains(p)))
                {
                    Point average = point.Points.Average();
                    //if (average != null)
                    {
                        point.MoveTo(average, false);
                    }
                });
        }

        private void CenterPoint(Point point, double precision, bool safe, bool edges = true)
        {
            if (point.Forced)
            {
                return;
            }

            int boundaryCount = point.BoundaryCount;
            if (boundaryCount == 0)
            {

                Point average = point.Points.Average();
                point.MoveTo(average, safe);
            }
            else if (boundaryCount == 1)
            {
                int ineqNumber = point.BoundaryFirstIndex;
                Point average = point.Points.Where(p => p.Boundary[ineqNumber]).Average();
                if (average != null)
                {
                    point.MoveTo(average, safe);
                    ProjectToSurface(point, precision, ineqNumber, safe);
                }
            }
            else if (edges && boundaryCount == 2)
            {
                int ineqNumber1 = point.BoundaryFirstIndex;
                int ineqNumber2 = point.BoundarySecondIndex;
                Point average = point.Points.Where(p => p.Boundary[ineqNumber1] && p.Boundary[ineqNumber2]).Average();
                if (average != null)
                {
                    point.MoveTo(average, safe);
                    ProjectToEdge(point, ineqNumber1, ineqNumber2, safe);
                }
            }
            else if (edges && boundaryCount == 3)
            {
                int ineqNumber1 = point.BoundaryFirstIndex;
                int ineqNumber2 = point.BoundarySecondIndex;
                int ineqNumber3 = point.BoundaryThirdIndex;
                ProjectToCorner(point, ineqNumber1, ineqNumber2, ineqNumber3, safe);
            }
        }

        private void Eval(Point p, int ineqNumber)
        {
            if (p.Values.ContainsKey(ineqNumber))
            {
                p.U = p.Values[ineqNumber];
            }
            else
            {
                p.U = Eval(p.X, p.Y, p.Z, ineqNumber);
                p.Values[ineqNumber] = p.U;
            }

            /*if (Math.Abs(p.U) < 1e-10)
            {
                p.U = 0;
            }*/
            }

        private double Eval(double x, double y, double z, int ineqNumber)
        {
            return ineqTreeBoxed.ExpressionList[ineqNumber](x, y, z);
        }

        public Func<double, double, double, double> EvalTotalFunc { get; set; }

        public double EvalTotal(Point p)
        {
            if (EvalTotalFunc != null)
                return EvalTotalFunc(p.X, p.Y, p.Z);

            double res = 0;
            for (int i = 0; i < ineqTreeBoxed.ExpressionList.Count; i++)
            {
                res += Math.Pow(Eval(p.X, p.Y, p.Z, i), 2);
            }
            return res;
        }

        public override Point AddPoint(double x, double y, double z)
        {
            Point p = new Point(x, y, z, ineqTreeBoxed.ExpressionList.Count);
            points.Add(p);
            return p;
        }

        private Point CreatePoint(double x, double y, double z)
        {
            Point p = new Point(x, y, z, ineqTreeBoxed.ExpressionList.Count);
            return p;
        }

        public Tetrahedron CreateTetrahedron(Point p0, Point p1, Point p2, Point p3)
        {
            Tetrahedron t = new Tetrahedron(p0, p1, p2, p3, boundaryCount, domainCount);
            return t;
        }


        private void ResolveMeshApriori(IneqTree.IneqNode node, int domaiNumber)
        {
            if (domaiNumber == 0)
                currentDomaiNumber = 0;

            if (node.NodeType == IneqTree.NodeType.NodeExpression)
            {
                ResolveIneqApriori(node.ExpressionIndex, domaiNumber);
            }
            else
            {
                int domaiNumberLeft = ++currentDomaiNumber;
                int domaiNumberRight = ++currentDomaiNumber;

                ResolveMeshApriori(node.Left, domaiNumberLeft);
                ResolveMeshApriori(node.Right, domaiNumberRight);

                int[] expressionIndexes = node.ExpressionIndexes.ToArray();

                Parallel.ForEach(Tetrahedrons, t =>
                {
                    if (node.NodeType == IneqTree.NodeType.NodeOr)
                    {
                        if (t.IsInDomain(domaiNumberLeft) || t.IsInDomain(domaiNumberRight))
                        {
                            t.IsIn[domaiNumber] = true;
                            t.IsOnBoundary[domaiNumber] = false;
                        }
                        else if (t.IsOnBoundaryDomain(domaiNumberLeft) || t.IsOnBoundaryDomain(domaiNumberRight))
                        {
                            t.IsIn[domaiNumber] = false;
                            t.IsOnBoundary[domaiNumber] = true;
                        }
                        else
                        {
                            t.IsIn[domaiNumber] = false;
                            t.IsOnBoundary[domaiNumber] = false;
                        }
                    }
                    if (node.NodeType == IneqTree.NodeType.NodeAnd)
                    {
                        if (t.IsOutDomain(domaiNumberLeft) || t.IsOutDomain(domaiNumberRight))
                        {
                            t.IsIn[domaiNumber] = false;
                            t.IsOnBoundary[domaiNumber] = false;
                        }
                        else if (t.IsOnBoundaryDomain(domaiNumberLeft) || t.IsOnBoundaryDomain(domaiNumberRight))
                        {
                            t.IsIn[domaiNumber] = false;
                            t.IsOnBoundary[domaiNumber] = true;
                        }
                        else
                        {
                            t.IsIn[domaiNumber] = true;
                            t.IsOnBoundary[domaiNumber] = false;
                        }
                    }
                    if (!t.IsOnBoundary[domaiNumber])
                    {
                        foreach (int i in expressionIndexes)
                        {
                            t.Boundary[i] = false;
                        }
                    }
                });
            }
        }

        private void ResolveIneqApriori(int ineqNumber, int domaiNumber)
        {
            if (OnProgress != null)
                OnProgress(0.5 * (double)(ineqNumber + 1) / (ineqTreeBoxed.ExpressionList.Count));

            Parallel.ForEach(Points, p =>
            {
                Eval(p, ineqNumber);
            });

            Parallel.ForEach(Tetrahedrons, t =>
            {
                if (t.Points.All(p => p.U < 0))
                {
                    t.IsIn[domaiNumber] = true;
                    t.IsOnBoundary[domaiNumber] = false;
                    t.Boundary[ineqNumber] = false;
                }
                else if (t.Points.All(p => p.U > 0))
                {
                    t.IsIn[domaiNumber] = false;
                    t.IsOnBoundary[domaiNumber] = false;
                    t.Boundary[ineqNumber] = false;
                }
                else
                {
                    t.IsIn[domaiNumber] = false;
                    t.IsOnBoundary[domaiNumber] = true;
                    t.Boundary[ineqNumber] = true;
                }
            });
        }

        private void SpreadTetrahedronsBoundaryFlags()
        {
            foreach (Tetrahedron t in Tetrahedrons.AsParallel().Where(t => t.BoundaryCount > 0))
            {
                foreach (Point p in t.Points)
                {
                    for (int i = 0; i < p.Boundary.Length; i++)
                    {
                        p.Boundary[i] = p.Boundary[i] || t.Boundary[i];
                    }
                }
            }

            foreach (Tetrahedron t in Tetrahedrons.AsParallel())
            {
                for (int i = 0; i < t.Boundary.Length; i++)
                {
                    t.Boundary[i] = t.Boundary[i] || t.Points.Any(p => p.Boundary[i]);
                }
            }

            Points.AsParallel().ForAll(p => p.Boundary.SetAll(false));
        }

        public void RefineBoundaryTriangle(Triangle t0)
        {
            LinkedList<Edge> refineEdges = new LinkedList<Edge>();
            Nullable<Triangle> t = t0;

            Edge e = t0.Edges().OrderByDescending(ee => ee.SqrLength).First();

            while (t != null)
            {
                refineEdges.AddFirst(e);

                Triangle[] trians = e
                    .P1.Tetrahedrons.Intersect(e.P2.Tetrahedrons)
                    .SelectMany(tt => tt.Triangles())
                    .Where(tr => tr.Boundary)
                    .Where(tr => tr.Contains(e.P1) && tr.Contains(e.P2) && !tr.Equals(t.Value))
                    /*.GroupBy(tr => tr)
                    .Where(gr => gr.Count() == 1)
                    .Select(gr => gr.Single())*/
                    .ToArray();

                if (trians.Length == 1)
                {
                    Edge e1 = trians[0].Edges().OrderByDescending(ee => ee.SqrLength).First();
                    if (!e1.Equals(e))
                    {
                        e = e1;
                        t = trians[0];
                    }
                    else
                        t = null;
                }
                else
                    t = null;
            }

            foreach (Edge ee in refineEdges)
            {
                var p = DivideEdge(ee, -1, (ee.P1 + ee.P2) / 2);

                if (p.BoundaryCount == 1)
                {
                    ProjectToSurface(p, 100, p.BoundaryFirstIndex, true);
                }
                else if (p.BoundaryCount == 2)
                {
                    ProjectToEdge(p, p.BoundaryFirstIndex, p.BoundarySecondIndex, true);
                }
            }
        }

        public void RefineBoundaryTriangles(IEnumerable<Triangle> triangleList)
        {
            foreach (Triangle t in triangleList.ToArray())
            {
                if (t.Valid)
                    RefineBoundaryTriangle(t);
            }
        }

        public int CheckQuality(double minQuality, bool allowDivide, Tetrahedron[] badTetra = null)
        {
            if (badTetra == null)
            {
                badTetra = Tetrahedrons
                    .AsParallel()
                    .Where(t => t.Quality < minQuality)
                    .OrderBy(t => t.Quality)
                    .ToArray();
            }

            foreach (var t in badTetra)
            {
                if (!Tetrahedrons.Contains(t) || t.Quality >= minQuality)
                    continue;

                var divideAndCollapse = t.Edges().Select(e =>
                        {
                            Point toPoint = null;
                            double minPQuality = 0;

                            TryDivideAndCollapse(e, ref toPoint, ref minPQuality, allowDivide);

                            return new { Edge = e, MinQuality = minPQuality, ToPoint = toPoint };
                        }
                    )
                    .OrderByDescending(ee => ee.MinQuality)
                    .First();

                if (divideAndCollapse.MinQuality > 0)
                {
                    Point midPoint = DivideEdge(divideAndCollapse.Edge, -1, (divideAndCollapse.Edge.P1 + divideAndCollapse.Edge.P2) / 2);
                    if (divideAndCollapse.ToPoint != null)
                    {
                        CollapseEdge(new Edge(midPoint, divideAndCollapse.ToPoint), divideAndCollapse.ToPoint, false); //partial Jiggle is problematic
                    }
                    else
                    {
                        //Jiggle(2, midPoint.Points.SelectMany(p => p.Points)); //partial Jiggle is problematic
                    }
                }
                else
                {
                    Edge shortEdge = t.Edges().Where(e => e.P1.HasAllBoundary(e.P2) || e.P2.HasAllBoundary(e.P1)).OrderBy(ee => ee.SqrLength).FirstOrDefault();

                    if (shortEdge.P1 != null && Math.Sqrt(shortEdge.SqrLength) < D / 2.0)
                    {
                        CollapseEdge(shortEdge, null, false); //partial Jiggle is problematic
                    }
                }
            }

            return badTetra.Length;
        }


        private void TryDivideAndCollapse(Edge e, ref Point toPoint, ref double minQuality, bool allowDivide)
        {
            Point midPoint = (e.P1 + e.P2) / 2;

            HashSet<Point> fanPoints = new HashSet<Point>();
            HashSet<Tetrahedron> fanTetrahedrons = new HashSet<Tetrahedron>();

            double origMinQuality = 1000;
            double pMinQuality = 1000;

            for (int i = 0; i < midPoint.Boundary.Length; i++)
            {
                midPoint.Boundary[i] = e.P1.Boundary[i] && e.P2.Boundary[i];
            }

            foreach (Tetrahedron t in e.P1.Tetrahedrons.Intersect(e.P2.Tetrahedrons))
            {
                Point[] p = t.Points.Where(pp => pp != e.P1 && pp != e.P2).ToArray();

                fanPoints.Add(p[0]);
                fanPoints.Add(p[1]);

                Tetrahedron tt = new Tetrahedron(p[0], p[1], e.P1, midPoint, boundaryCount, domainCount, true);
                fanTetrahedrons.Add(tt);
                tt = new Tetrahedron(p[0], p[1], e.P2, midPoint, boundaryCount, domainCount, true);
                fanTetrahedrons.Add(tt);

                if (t.Quality < origMinQuality)
                {
                    origMinQuality = t.Quality;
                }
            }

            toPoint = null;
            minQuality = 0;

            foreach (Point p in fanPoints)
            {
                if (!midPoint.HasAllBoundary(p))
                    continue;

                pMinQuality = 1000;

                midPoint.MoveTo(p, false);
                foreach (var t in fanTetrahedrons.Where(t => !t.Points.Contains(p)))
                {
                    if (!t.CheckVolume(0.15d))
                    {
                        pMinQuality = -1;
                    }
                    else if (t.Quality < pMinQuality)
                    {
                        pMinQuality = t.Quality;
                    }
                }

                if (pMinQuality > minQuality)
                {
                    minQuality = pMinQuality;
                    toPoint = p;
                }
            }

            if (allowDivide && minQuality <= origMinQuality)
            {
                Point average = fanPoints.Union(e).Average();
                midPoint.MoveTo(average, false);
                pMinQuality = 1000;
                foreach (var t in fanTetrahedrons)
                {
                    if (!t.CheckVolume(0.05d))
                    {
                        pMinQuality = -1;
                    }
                    else if (t.Quality < pMinQuality)
                    {
                        pMinQuality = t.Quality;
                    }
                }

                if (pMinQuality > minQuality)
                {
                    minQuality = pMinQuality;
                    toPoint = null;
                }
            }

            if (minQuality < origMinQuality)
            {
                minQuality = 0;
                toPoint = null;
            }
        }


        public int CheckBoundaryQuality(double minQuality, bool allowDivide)
        {
            var badTriangles = Tetrahedrons.SelectMany(t => t.Triangles().Where(tr => tr.Boundary))
                                .Where(t => t.Quality < minQuality)
                                .ToArray();


            foreach (var t in badTriangles)
            {
                if (!t.Valid || t.Quality >= minQuality)
                    continue;

                var divideAndCollapse = t.Edges().Select(e =>
                        {
                            Point toPoint = null;
                            double minPQuality = 0;

                            TryDivideAndCollapseBoundary(e, ref toPoint, ref minPQuality, allowDivide);

                            return new { Edge = e, MinQuality = minPQuality, ToPoint = toPoint };
                        }
                    )
                    .OrderByDescending(ee => ee.MinQuality)
                    .First();

                if (divideAndCollapse.MinQuality > 0)
                {
                    Point midPoint = DivideEdge(divideAndCollapse.Edge, -1, (divideAndCollapse.Edge.P1 + divideAndCollapse.Edge.P2) / 2);
                    if (divideAndCollapse.ToPoint != null)
                    {
                        CollapseEdge(new Edge(midPoint, divideAndCollapse.ToPoint), divideAndCollapse.ToPoint, false);
                    }
                }
                //else
                //{
                //    Edge shortEdge = t.Edges().Where(e => e.P1.HasAllBoundary(e.P2) || e.P2.HasAllBoundary(e.P1)).OrderBy(ee => ee.SqrLength).FirstOrDefault();

                //    if (shortEdge.P1 != null && Math.Sqrt(shortEdge.SqrLength) < D / 2.0)
                //    {
                //        CollapseEdge(shortEdge, null, false); //partial Jiggle is problematic
                //    }
                //}
            }

            return badTriangles.Length;
        }


        private void TryDivideAndCollapseBoundary(Edge e, ref Point toPoint, ref double minQuality, bool allowDivide)
        {
            var points = e.P1.Tetrahedrons.Intersect(e.P2.Tetrahedrons)
                .SelectMany(t => t.Points.Where(pp => pp != e.P1 && pp != e.P2))
                .GroupBy(p => p)
                .Where(g => g.Count() == 1)
                .Select(g => g.Key)
                .ToArray();


            if (points.Length != 2)
            {
                minQuality = 0;
                toPoint = null;
                return;
            }

            Point midPoint = (e.P1 + e.P2) / 2;
            for (int i = 0; i < midPoint.Boundary.Length; i++)
            {
                midPoint.Boundary[i] = e.P1.Boundary[i] && e.P2.Boundary[i];
            }

            Triangle t1 = new Triangle(e.P1, e.P2, points[0]);
            Triangle t2 = new Triangle(e.P1, e.P2, points[1]);

            Triangle nt1 = new Triangle(e.P1, points[0], points[1]);
            Triangle nt2 = new Triangle(e.P2, points[0], points[1]);


            double origMinQuality = Math.Min(t1.Quality, t2.Quality);
            minQuality = Math.Min(nt1.Quality, nt2.Quality);


            if (minQuality > origMinQuality)
            {
                if (midPoint.HasAllBoundary(points[0]))
                {
                    toPoint = points[0];
                }
                else if (midPoint.HasAllBoundary(points[1]))
                {
                    toPoint = points[1];
                }
                else
                {
                    minQuality = 0;
                    toPoint = null;
                    return;
                }

                HashSet<Tetrahedron> fanTetrahedrons = new HashSet<Tetrahedron>();
                foreach (Tetrahedron t in e.P1.Tetrahedrons.Intersect(e.P2.Tetrahedrons))
                {
                    Point[] p = t.Points.Where(pp => pp != e.P1 && pp != e.P2).ToArray();

                    Tetrahedron tt = new Tetrahedron(p[0], p[1], e.P1, midPoint, boundaryCount, domainCount, true);
                    fanTetrahedrons.Add(tt);
                    tt = new Tetrahedron(p[0], p[1], e.P2, midPoint, boundaryCount, domainCount, true);
                    fanTetrahedrons.Add(tt);
                }

                midPoint.MoveTo(toPoint, false);
                var tp = toPoint;
                foreach (var t in fanTetrahedrons.Where(t => !t.Points.Contains(tp)))
                {
                    if (!t.CheckVolume(0.01d))
                    {
                        minQuality = 0;
                        toPoint = null;
                    }
                }
            }
            else
            {
                minQuality = 0;
                toPoint = null;
            }

            /*
                        toPoint = null;
                        minQuality = 0;

                        foreach (Point p in fanPoints)
                        {
                            if (!midPoint.HasAllBoundary(p))
                                continue;

                            pMinQuality = 1000;

                            midPoint.MoveTo(p, false);
                            foreach (var t in fanTetrahedrons.Where(t => !t.Points.Contains(p)))
                            {
                                if (!t.CheckVolume(0.15d))
                                {
                                    pMinQuality = -1;
                                }
                                else if (t.Quality < pMinQuality)
                                {
                                    pMinQuality = t.Quality;
                                }
                            }

                            if (pMinQuality > minQuality)
                            {
                                minQuality = pMinQuality;
                                toPoint = p;
                            }
                        }

                        if (allowDivide && minQuality <= origMinQuality)
                        {
                            Point average = fanPoints.Union(e).Average();
                            midPoint.MoveTo(average, false);
                            pMinQuality = 1000;
                            foreach (var t in fanTetrahedrons)
                            {
                                if (!t.CheckVolume(0.05d))
                                {
                                    pMinQuality = -1;
                                }
                                else if (t.Quality < pMinQuality)
                                {
                                    pMinQuality = t.Quality;
                                }
                            }

                            if (pMinQuality > minQuality)
                            {
                                minQuality = pMinQuality;
                                toPoint = null;
                            }
                        }

                        if (minQuality < origMinQuality)
                        {
                            minQuality = 0;
                            toPoint = null;
                        }
                        */
        }


        public void CheckCurvatureQuality(int count = 1)
        {
            var good = count == 1 ? null : new HashSet<Triangle>();

            for (int j = 0; j < count; j++)
            {
                List<Triangle> refList = new List<Triangle>();

                var edges = Edges
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

                    ProjectToEdge(ee.p, ee.bf1, ee.bf2, false);

                    double dist = origp.Distance(ee.p);

                    if (dist >= ee.length / 50.0d)
                    {
                        var e = ee.e;

                        var trians = e
                            .P1.Tetrahedrons.Intersect(e.P2.Tetrahedrons)
                            .SelectMany(tt => tt.Triangles())
                            .Where(tr => tr.Boundary)
                            .Where(tr => tr.Contains(e.P1) && tr.Contains(e.P2));
                        /*.GroupBy(tr => tr)
                        .Where(gr => gr.Count() == 1)
                        .Select(gr => gr.Single());*/

                        refList.AddRange(trians);

                    }
                }

                var centerPoints = Tetrahedrons.AsParallel().SelectMany(t => t.Triangles()
                                .Where(tr => tr.BoundaryCount == 1 && tr.Boundary && (count == 1 || !good.Contains(tr))))
                                .Select(tr => new
                                {
                                    bf = tr.CommonBoundaryFlag.Value,
                                    p = tr.Average(),
                                    tr,
                                    maxLength = Math.Max(Math.Max(tr.P1.Distance(tr.P2), tr.P1.Distance(tr.P3)), tr.P2.Distance(tr.P3))
                                });

                //foreach (var cp in centerPoints)
                centerPoints.ForAll(cp =>
                {
                    Point origp = new Point(cp.p.X, cp.p.Y, cp.p.Z);

                    ProjectToSurface(cp.p, 100, cp.bf, false);

                    double dist = origp.Distance(cp.p);

                    if (dist >= cp.maxLength / 50.0d)
                    {
                        lock (refList)
                        {
                            refList.Add(cp.tr);
                        }
                    }
                    else if (count > 1)
                    {
                        lock (good)
                        {
                            good.Add(cp.tr);
                        }
                    }
                });

                RefineBoundaryTriangles(refList);
            }

            DeleteLonelyPoints();

            Jiggle(3, false);
        }


        /*private void PrecalsCorners()
        {
            List<Point> corners = new List<Point>(20);

            foreach (Tetrahedron t in Tetrahedrons.AsParallel().Where(t => t.BoundaryCount >= 3))
            {
                int[] ineqNumbers = t.Boundary.Cast<bool>().Select((b, i) => new { IsSet = b, index = i }).Where(bf => bf.IsSet).Take(3).Select(bi => bi.index).ToArray();

                Point p = t.Points.Average();

                if (ProjectToCorner(p, ineqNumbers[0], ineqNumbers[1], ineqNumbers[2], false) && !corners.Any(cp => cp.Distance(p) < D /10))
                {
                    Point p1 = t.Points.OrderBy(pp => pp.Distance(p)).First();

                    if (p1.MoveTo(p, true))
                    {
                        corners.Add(p);

                        p1.Boundary[ineqNumbers[0]] = true;
                        p1.Boundary[ineqNumbers[1]] = true;
                        p1.Boundary[ineqNumbers[2]] = true;

                        foreach(var pp in p1.Points)
                        {
                            CenterPoint(pp, 100, true);
                        }
                    }
                } 
            }               
        }*/

        public void CheckTopology()
        {
            List<Tetrahedron> tetras = null;

            do
            {
                tetras = Tetrahedrons.Where(t => t.P0.BoundaryCount > 0 && t.CommonBoundaryCount > 0).ToList();
                foreach (var tetra in tetras)
                {
                    DeleteTetrahedron(tetra);
                }
            }
            while (tetras.Count > 0);

            foreach (var p in Points.Where(pp => pp.BoundaryCount == 2))
            {
                var b1 = p.BoundaryFirstIndex;
                var b2 = p.BoundarySecondIndex;

                if (p.Points.Where(pp => pp.BoundaryCount > 0).All(pp => pp.Boundary[b1]))
                {
                    p.Boundary[b2] = false;
                }

                if (p.Points.Where(pp => pp.BoundaryCount > 0).All(pp => pp.Boundary[b2]))
                {
                    p.Boundary[b1] = false;
                }
            }

            foreach (var p in Points.Where(pp => pp.BoundaryCount > 2))
            {
                foreach (int ineqNumber in p.Boundary.Cast<bool>().Select((b, i) => new { b = b, i = i }).Where(bi => bi.b).Select(bi => bi.i).ToArray())
                {
                    if (!p.Points.Any(p1 => p1.Boundary[ineqNumber] && p1.BoundaryCount < p.BoundaryCount))
                        p.Boundary[ineqNumber] = false;
                }
            }

            DeleteLonelyPoints();

            return;
        }

        //*******************************************************************************
        /*public Point ForcePoint(Point p, bool fixedPoint, Point p1)
        {
            if (p.Forced)
            {
                return p;
            }

            var np = Points.OrderBy(pp => pp.Distance(p)).First();

            if (np.Forced && np.Distance(p) <= 1e-4)
            {
                return np;
            }

            if (!np.Forced && np.Movable && np.MoveTo(p, true))
            {
                np.Forced = true;
                np.Movable = false;

                return np;
            }

            var t = Tetrahedrons.FirstOrDefault(tt => tt.ContainsPoint(p));

            if (fixedPoint)
            {
                p = AddPoint(p.X, p.Y, p.Z);
                AddTetrahedron(t.P0, t.P1, t.P2, p);
                AddTetrahedron(t.P0, t.P1, p, t.P3);
                AddTetrahedron(t.P0, p, t.P2, t.P3);
                AddTetrahedron(p, t.P1, t.P2, t.P3);

                DeleteTetrahedron(t);

                p.Forced = true;
                p.Movable = false;

            }
            else
            {
                var intRes = t.Triangles().Select(trian => new { trian, res = trian.IntersectLine(p, p1)}).Where(r => r.res.Point != null && r.res.Point != p1 && !r.res.Point.Forced).OrderBy(r => p.Distance(r.res.Point)).FirstOrDefault();

                if (intRes != null && intRes.res.Type == "inside")
                {
                    var tr = intRes.trian;
                    var intP = intRes.res.Point;

                    var t1 = tr.P1.Tetrahedrons.Intersect(tr.P2.Tetrahedrons).Intersect(tr.P3.Tetrahedrons).Where(tt => tt != t).Single();

                    p = AddPoint(intP.X, intP.Y, intP.Z);

                    var p4 = t.Points.Single(pp => pp != tr.P1 && pp != tr.P2 && pp != tr.P3);

                    AddTetrahedron(tr.P1, tr.P2, p, p4);
                    AddTetrahedron(tr.P2, tr.P3, p, p4);
                    AddTetrahedron(tr.P3, tr.P1, p, p4);

                    DeleteTetrahedron(t);

                    p4 = t1.Points.Single(pp => pp != tr.P1 && pp != tr.P2 && pp != tr.P3);

                    AddTetrahedron(tr.P1, tr.P2, p, p4);
                    AddTetrahedron(tr.P2, tr.P3, p, p4);
                    AddTetrahedron(tr.P3, tr.P1, p, p4);

                    DeleteTetrahedron(t1);


                    p.Forced = true;
                    p.Movable = false;
                }
                else if (intRes != null && intRes.res.Type == "vertex")
                {
                    p = intRes.res.Point;
                    p.Forced = true;
                    p.Movable = false;
                }
                else if (intRes != null && intRes.res.Type == "edge")
                {
                    var intP = intRes.res.Point;
                    p = new Point(intP.X, intP.Y, intP.Z, ineqTreeBoxed.ExpressionList.Count);

                    DivideEdge(intRes.res.Edge.Value, -1, p);

                    p.Forced = true;
                    p.Movable = false;
                }
                else
                {
                    bool bad = true;
                }
            }

            return p;
        }*/

        /*public void ForceEdge(Point p1, Point p2, IneqTree ineq)
        {
            ForceEdge(ref p1, ref p2, ineq, 0, true);
        }

        private void ForceEdge(ref Point p1, ref Point p2, IneqTree ineq, int l, bool fixedPoints)
        {
            if (l > 8)
            {
                return;
            }

            p1 = ForcePoint(p1, fixedPoints, p2);
            p2 = ForcePoint(p2, fixedPoints, p1);

            if (!p1.Forced || !p2.Forced)
            {
                return;
            }

            if(!p1.Points.Contains(p2) && p1.Forced && p2.Forced && (p1 != p2))
            {
                var cp = (p1 + p2) / 2;

                ForceEdge(ref p1, ref cp, ineq, l + 1, false);
                ForceEdge(ref p2, ref cp, ineq, l + 1, false);
            }
            else if (p1.Forced && p2.Forced && (p1 != p2) && ineq != null && ineq.Root.NodeType == IneqTree.NodeType.NodeExpression)
            {
                int ineqNumber = ineqTreeBoxed.ExpressionList.IndexOf(ineq.Root.Expression);

                var np1 = p1;
                var np2 = p2;
                var points = p1.Tetrahedrons.Intersect(p2.Tetrahedrons).SelectMany(t => t.Points).Where(p => p != np1 && p != np2).Distinct();

                foreach (var p in points)
                {
                    Eval(p, ineqNumber);
                }

                if (points.All(p => p.U >= 0))
                {
                    var edges = p1.Tetrahedrons.Intersect(p2.Tetrahedrons).SelectMany(t => t.Edges()).Where(e => (e.P1 != np1 && e.P1 != np2) || (e.P2 != np1 && e.P2 != np2)).Distinct();

                    foreach (var e in edges.ToArray())
                    {
                        if (!e.Valid)
                            continue;

                        Eval(e.P1, ineqNumber);
                        Eval(e.P2, ineqNumber);

                        if (e.P1.U * e.P2.U < 0)
                            continue;

                        int count = 10;
                        var f = Enumerable.Range(1, count - 1).Select(i =>
                        {
                            var mp = e.P1 + (((double)i) / count) * (e.P2 - e.P1);
                            Eval(mp, ineqNumber);
                            return new { i, mp };
                        }).Where(it => it.mp.U < 0).OrderBy(it => Math.Abs(it.i - count / 2)).FirstOrDefault();

                        if (f != null)
                        {
                            DivideEdge(e, -1, f.mp);
                            f.mp.Movable = false;
                            foreach (var pp in f.mp.Points)
                            {
                                pp.Movable = false;
                            }
                        }
                    }
                }
            }

            foreach (var pp in p1.Points)
            {
                pp.Movable = false;
            }

            foreach (var pp in p2.Points)
            {
                pp.Movable = false;
            }
        }
        */


        //********************************


        //private HashSet<Point> forcedPoints = new HashSet<Point>();
        //public void ForceEdge(Point p1, Point p2, double[] n1, double[] n2)
        //{
        //    forcedPoints.Clear();
        //    ForceEdge(ref p1, ref p2, n1, n2, 0, true);

        //    var edges = forcedPoints.SelectMany(p => p.Tetrahedrons).SelectMany(t => t.Edges()).Distinct();

        //    foreach (var e in edges.ToArray())
        //    {
        //        if (!e.Valid)
        //            continue;

        //        /*TestEdgePlane(e, p1, n1);
        //        TestEdgePlane(e, p1, n2);*/

        //    }

        //    foreach (var p in forcedPoints.SelectMany(pp => pp.Points))
        //    { 
        //        p.Movable = false;
        //    }
        //}

        //private void ForceEdge(ref Point p1, ref Point p2, double[] n1, double[] n2, int l, bool fixedPoints)
        //{
        //    if (l > 8)
        //    {
        //        return;
        //    }

        //    p1 = ForcePoint(p1, fixedPoints, p2);
        //    p2 = ForcePoint(p2, fixedPoints, p1);

        //    if (!p1.Forced || !p2.Forced)
        //    {
        //        return;
        //    }

        //    forcedPoints.Add(p1);
        //    forcedPoints.Add(p2);


        //    if (!p1.Points.Contains(p2) && (p1 != p2))
        //    {
        //        var cp = (p1 + p2) / 2;

        //        ForceEdge(ref p1, ref cp, n1, n2, l + 1, false);
        //        ForceEdge(ref p2, ref cp, n1, n2, l + 1, false);
        //    }
        //    /*else if (p1 != p2)
        //    {
        //        var np1 = p1;
        //        var np2 = p2;
        //        var edges = p1.Tetrahedrons.Intersect(p2.Tetrahedrons).SelectMany(t => t.Edges()).Where(e => (e.P1 != np1 && e.P1 != np2) || (e.P2 != np1 && e.P2 != np2)).Distinct();

        //        foreach (var e in edges.ToArray())
        //        {
        //            if (!e.Valid)
        //                continue;

        //            TestEdgePlane(e, p1, n1);
        //            TestEdgePlane(e, p1, n2);

        //        }
        //    }*/

        //    foreach (var pp in p1.Points)
        //    {
        //        pp.Movable = false;
        //    }

        //    foreach (var pp in p2.Points)
        //    {
        //        pp.Movable = false;
        //    }
        //}

        /*private void TestEdgePlane(Edge e, Point p, double[] n)
        {
            if(!e.Valid) 
                return;

            Func<double, double> f = t =>
            {
                var pt = (e.P1 + t * (e.P2 - e.P1)) - p;

                var res = pt.X * n[0] + pt.Y * n[1] + pt.Z * n[2];

                return res;
            };

            var f0 = f(0);
            var f1 = f(1);

            if (f0 * f1 < 0)
            {
                var t0 = f0 / (f0 - f1);
                var p0 = e.P1 + t0 * (e.P2 - e.P1);
                DivideEdge(e, -1, p0);
                p0.Movable = false;
            }
        }*/

        //************************* 

        public void ForceEdge(Point p1, Point p2, double[] n1, double[] n2)
        {
            /*p1 = ForcePoint(p1);
            p2 = ForcePoint(p2);*/

            var f1 = GetLinearFunc(p1, n1);
            var f2 = GetLinearFunc(p1, n2);
            var cutTetras = new HashSet<Tetrahedron>();

            var addTetras = Tetrahedrons.Where(tt => tt.ContainsPoint(p1)).ToArray();

            while (addTetras.Length > 0)
            {
                cutTetras.UnionWith(addTetras);

                addTetras = addTetras
                    .SelectMany(t => t.Points)
                    .SelectMany(p => p.Tetrahedrons)
                    .Where(t => !cutTetras.Contains(t) && IntersectTetra(t, f1) && IntersectTetra(t, f2) && t.Points.Any(p => BetweenPoints(p, p1, p2)))
                    .Distinct()
                    .ToArray();
            }

            foreach (var e in cutTetras.SelectMany(t => t.Edges()).Distinct().ToArray())
            {
                CutEdgePlane(e, p1, n1, cutTetras);
            }

            foreach (var e in cutTetras.SelectMany(t => t.Edges()).Distinct().ToArray())
            {
                CutEdgePlane(e, p1, n2, cutTetras);
            }
        }

        private FuncXYZ GetLinearFunc(Point p, double[] n)
        {
            return (x, y, z) =>
            {
                return (x - p.X) * n[0] + (y - p.Y) * n[1] + (z - p.Z) * n[2];
            };
        }

        private bool IntersectTetra(Tetrahedron t, FuncXYZ f)
        {
            var vals = t.Points.Select(p => f(p.X, p.Y, p.Z))/*.Where(v => v != 0)*/.ToArray();

            return !(vals.All(v => v > 0) || vals.All(v => v < 0));
        }

        private bool BetweenPoints(Point p, Point p1, Point p2)
        {
            // Vector p1 → p2
            double vx = p2.X - p1.X;
            double vy = p2.Y - p1.Y;
            double vz = p2.Z - p1.Z;

            // Vector p1 → p
            double wx = p.X - p1.X;
            double wy = p.Y - p1.Y;
            double wz = p.Z - p1.Z;

            // Dot products
            double dot = vx * wx + vy * wy + vz * wz;
            if (dot < 0) return false;                   // projection is before p1

            double len2 = vx * vx + vy * vy + vz * vz;
            if (dot > len2) return false;                // projection is after p2

            return true;                                 // projection lies between p1 and p2
        }

        private void CutEdgePlane(Edge e, Point p, double[] n, HashSet<Tetrahedron> cutTetras)
        {
            if (!e.Valid)
                return;

            Func<double, double> f = t =>
            {
                var pt = (e.P1 + t * (e.P2 - e.P1)) - p;

                var res = pt.X * n[0] + pt.Y * n[1] + pt.Z * n[2];

                return res;
            };

            var f0 = f(0);
            var f1 = f(1);

            if (Math.Abs(f0) < 1E-6)
            {
                e.P1.Movable = false;
            }
            else if (Math.Abs(f1) < 1E-6)
            {
                e.P2.Movable = false;
            }
            else if (f0 * f1 < 0)
            {
                var t0 = f0 / (f0 - f1);

                if (t0 < 0 || t0 > 1)
                {
                    return;
                }

                var p0 = e.P1 + t0 * (e.P2 - e.P1);

                if (t0 < 0.45 && e.P1.Movable && e.P1.MoveTo(p0, true))
                {
                    e.P1.Movable = false;
                }
                else if (t0 > 0.55 && e.P2.Movable && e.P2.MoveTo(p0, true))
                {
                    e.P2.Movable = false;
                }
                else
                {
                    foreach (Tetrahedron t in e.P1.Tetrahedrons.Intersect(e.P2.Tetrahedrons).ToArray())
                    {
                        cutTetras.Remove(t);
                    }
                    DivideEdge(e, -1, p0);
                    p0.Movable = false;
                    if (!e.P1.Movable && !e.P2.Movable)
                    {
                        p0.Forced = true;
                    }
                    cutTetras.UnionWith(p0.Tetrahedrons);
                }
            }
        }

        public Point ForcePoint(Point p)
        {
            if (p.Forced)
            {
                return p;
            }

            var np = Points.OrderBy(pp => pp.Distance(p)).First();

            if (np.Forced && np.Distance(p) <= 1e-6)
            {
                return np;
            }

            if (!np.Forced && np.Movable && np.MoveTo(p, true))
            {
                np.Forced = true;
                np.Movable = false;

                return np;
            }

            var t = Tetrahedrons.FirstOrDefault(tt => tt.ContainsPoint(p));

            p = AddPoint(p.X, p.Y, p.Z);
            AddTetrahedron(t.P0, t.P1, t.P2, p);
            AddTetrahedron(t.P0, t.P1, p, t.P3);
            AddTetrahedron(t.P0, p, t.P2, t.P3);
            AddTetrahedron(p, t.P1, t.P2, t.P3);

            DeleteTetrahedron(t);

            p.Forced = true;
            p.Movable = false;

            return p;
        }
    }
}