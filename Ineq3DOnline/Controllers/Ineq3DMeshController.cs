using Ineq3DOnline.Models;
using MeshData;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Web;
using System.Web.Mvc;

namespace Ineq3DOnline.Controllers
{
    public class Ineq3DMeshController : Controller
    {
        public ActionResult Index(string Mesh = "")
        {
            MeshSamples samples = new MeshSamples();

            if (String.IsNullOrEmpty(Mesh))
            {
                Mesh = samples.Samples.Skip(new Random().Next(samples.Samples.Count())).First();
                //Mesh = "Ball";
            }

            ViewBag.CurrentMesh = Mesh;
            return View(samples);
        }

        public ActionResult GetMesh(string Mesh = "")
        {
            IneqMesh ineqMesh = null;
            string path = null;

            MeshSamples samples = new MeshSamples();

            string FileName = Mesh + ".ply";
            path = System.IO.Path.Combine(Server.MapPath("~/Samples"), FileName);
            if (System.IO.File.Exists(path))
            {
                return File(path, "application/octet-stream");
            }

            ineqMesh = samples[Mesh];
            ineqMesh.Create();
            CheckQuality(ineqMesh);
            ineqMesh.DeleteLonelyPoints();

            string ply = GetPLY(ineqMesh);

            System.IO.File.WriteAllText(path, ply);
            return File(path, "application/octet-stream");
        }


        [HttpPost, ValidateInput(false)]
        public ActionResult SetIneqMesh(IneqMeshViewModel ineqMeshViewModel)
        {
            if (ModelState.IsValid)
            {
                try
                {
                    ineqMeshViewModel.SetIneqMesh();
                    Session["IneqMeshViewModel"] = ineqMeshViewModel;
                }
                catch (Exception exc)
                {
                    return Json(new { success = false, message = exc.Message });
                }

                return Json(new { success = true });
            }

            return Json(new { success = false, message = "Invalid data" });
        }

        public ActionResult GetCustomMesh()
        {
            IneqMeshViewModel ineqMeshViewModel = (IneqMeshViewModel)Session["IneqMeshViewModel"];
            IneqMesh ineqMesh = null;

            if (ineqMeshViewModel == null)
            {
                return Content(null);
            }

            if (!String.IsNullOrEmpty(ineqMeshViewModel.PLY))
            {
                return Content(ineqMeshViewModel.PLY);
            }

            ineqMesh = ineqMeshViewModel.IneqMesh;
            ineqMesh.Create();

            if (ineqMeshViewModel.Quality)
            {
                CheckQuality(ineqMesh);
            }

            if (ineqMeshViewModel.CurvatureQuality)
            {
                CheckCurvatureQuality(ineqMesh);
            }

            string ply = GetPLY(ineqMesh);

            ineqMeshViewModel.PLY = ply;
            return Content(ply);
        }

        public ActionResult GetCustomMeshImproveQuality(bool boundary)
        {
            IneqMeshViewModel ineqMeshViewModel = (IneqMeshViewModel)Session["IneqMeshViewModel"];
            IneqMesh ineqMesh = null;

            if (ineqMeshViewModel == null || ineqMeshViewModel.IneqMesh == null || ineqMeshViewModel.IneqMesh.Tetrahedrons.Count == 0)
            {
                return Content(null);
            }

            ineqMesh = ineqMeshViewModel.IneqMesh;

            if (boundary)
            {
                CheckBoundaryQuality(ineqMesh);
            }
            else
            {
                CheckQuality(ineqMesh);
            }

            string ply = GetPLY(ineqMesh);

            ineqMeshViewModel.PLY = ply;
            return Content(ply);
        }

        public ActionResult GetCustomMeshImproveCurvatureQuality()
        {
            IneqMeshViewModel ineqMeshViewModel = (IneqMeshViewModel)Session["IneqMeshViewModel"];
            IneqMesh ineqMesh = null;

            if (ineqMeshViewModel == null || ineqMeshViewModel.IneqMesh == null || ineqMeshViewModel.IneqMesh.Tetrahedrons.Count == 0)
            {
                return Content(null);
            }

            ineqMesh = ineqMeshViewModel.IneqMesh;

            CheckCurvatureQuality(ineqMesh);
            CheckCurvatureQuality(ineqMesh);

            ineqMesh.DeleteLonelyPoints();

            string ply = GetPLY(ineqMesh);

            ineqMeshViewModel.PLY = ply;
            return Content(ply);
        }

        public ActionResult GetCustomMeshJiggle()
        {
            IneqMeshViewModel ineqMeshViewModel = (IneqMeshViewModel)Session["IneqMeshViewModel"];
            IneqMesh ineqMesh = null;

            if (ineqMeshViewModel == null || ineqMeshViewModel.IneqMesh == null || ineqMeshViewModel.IneqMesh.Tetrahedrons.Count == 0)
            {
                return Content(null);
            }

            ineqMesh = ineqMeshViewModel.IneqMesh;

            ineqMesh.Jiggle(3);

            ineqMesh.DeleteLonelyPoints();

            string ply = GetPLY(ineqMesh);

            ineqMeshViewModel.PLY = ply;
            return Content(ply);
        }

        public ActionResult GetSampleFormula(int sampleIndex)
        {
            var tmp = IneqMeshViewModel.DefaultModel(sampleIndex);

            return Json(new
            {
                formula = tmp.Formula
            });
        }

        public ActionResult GetSampleUFunc(int sampleUFuncIndex)
        {
            var tmp = IneqMeshViewModel.DefaultModel(-1, sampleUFuncIndex);

            return Json(new
            {
                ufunc = tmp.UFunc
            });
        }



        double minQuality = 0.25d;
        private void CheckQuality(IneqMesh ineqMesh)
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
        private void CheckBoundaryQuality(IneqMesh ineqMesh)
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

        private void CheckCurvatureQuality(IneqMesh ineqMesh)
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

        private string GetPLY(IneqMesh ineqMesh)
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
    }
}

//if (Mesh == "Cone 1")
//{
//    Point p = new Point(0.12, -0.171, 0.037);
//    Point np = null;
//    double d = Double.MaxValue;
//    foreach (var mp in ineqMesh.Points.Where(pp => pp.BoundaryCount == 1))
//    {
//        if (mp.Distance(p) < d)
//        {
//            d = mp.Distance(p);
//            np = mp;
//        }
//    }

//    np.MoveTo(p, false);
//    np.Boundary[np.Boundary.Length - 1] = true;
//    ineqMesh.Jiggle(3);

//    ineqMesh.RefineBoundaryTriangles(np.Tetrahedrons.SelectMany(t => t.Triangles()).Where(t => t.All(pp => pp.BoundaryCount > 0)));
//    ineqMesh.Jiggle(3);
//    ineqMesh.RefineBoundaryTriangles(ineqMesh.Triangles.Where(tt => tt.All(pp => pp.Boundary[0])));
//    ineqMesh.Jiggle(3);
//}

//ineqMesh.RefineBoundaryTriangles(ineqMesh.Triangles.Where(tt => tt.All(p => p.Boundary[5])));
//ineqMesh.Jiggle(3);

