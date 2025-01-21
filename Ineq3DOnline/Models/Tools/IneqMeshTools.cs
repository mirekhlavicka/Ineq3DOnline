using MeshData;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Web;

namespace Ineq3DOnline.Models
{
    public class IneqMeshTools
    {
        const double minQuality = 0.25d;
        public static void CheckQuality(IneqMesh ineqMesh)
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

        public static void CheckBoundaryQuality(IneqMesh ineqMesh)
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

        public static void CheckCurvatureQuality(IneqMesh ineqMesh)
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

            var centerPoints = ineqMesh.Tetrahedrons.SelectMany(t => t.Triangles()
                            .Where(tr => tr.BoundaryCount == 1 && tr.Boundary))
                            .Select(tr => new
                            {
                                bf = tr.CommonBoundaryFlag.Value,
                                p = tr.Average(),
                                tr,
                                maxLength = Math.Max(Math.Max(tr.P1.Distance(tr.P2), tr.P1.Distance(tr.P3)), tr.P2.Distance(tr.P3))
                            });

            foreach (var cp in centerPoints)
            {
                Point origp = new Point(cp.p.X, cp.p.Y, cp.p.Z);

                ineqMesh.ProjectToSurface(cp.p, 100, cp.bf, false);

                double dist = origp.Distance(cp.p);

                if (dist >= cp.maxLength / 50.0d)
                {
                    refList.Add(cp.tr);
                }
            }

            ineqMesh.RefineBoundaryTriangles(refList);

            ineqMesh.DeleteLonelyPoints();

            ineqMesh.Jiggle(3);
        }
    }
}