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
        string plyFormat =
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


        string[] colors = {
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
            IneqMeshViewModel ineqMeshViewModel = null;
            IneqMesh ineqMesh = null;
            string path = null;

            if (Mesh == "[custom]")
            {
                ineqMeshViewModel = (IneqMeshViewModel)Session["IneqMeshViewModel"];

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
                    ineqMesh.DeleteLonelyPoints();
                }
            }
            else
            {
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

            StringBuilder sbPoints = new StringBuilder();
            StringBuilder sbTriangles = new StringBuilder();
            int pointsCount = 0;

            var boundaryTriangles = ineqMesh.Tetrahedrons.SelectMany(t => t.Triangles().Where(tr => /*tr.BoundaryCount == 1 &&*/ tr.P1.Tetrahedrons.Intersect(tr.P2.Tetrahedrons).Intersect(tr.P3.Tetrahedrons).Count() == 1)).ToArray();

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


            string ply = String.Format(plyFormat, pointsCount, boundaryTriangles.Length, sbPoints, sbTriangles);

            if (path != null)
            {
                System.IO.File.WriteAllText(path, ply);
                return File(path, "application/octet-stream");
            }
            else if (ineqMeshViewModel != null)
            {
                ineqMeshViewModel.PLY = ply;
                return Content(ply);
            }
            else
            {
                return Content(null);
            }
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

        double minQuality = 0.3d;
        private void CheckQuality(IneqMesh mesh)
        {
            if (mesh == null)
                return;


            int c = 0;
            Dictionary<int, int> counts = new Dictionary<int, int>();
            do
            {
                c = mesh.CheckQuality(minQuality, false);

                if (!counts.Keys.Contains(c))
                    counts[c] = 1;
                else
                    counts[c]++;


            }
            while (c != 0 && counts[c] < 3);

            if (c != 0)
            {
                c = mesh.CheckQuality(minQuality, true);
            }

            return;
        }
    }
}
