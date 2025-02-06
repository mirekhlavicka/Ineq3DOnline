using Ineq3DOnline.Models;
using MeshData;
using Newtonsoft.Json.Serialization;
using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Web;
using System.Web.Mvc;
using WebGrease.Activities;

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

        public ActionResult SetSampleMesh(string Mesh = "")
        {
            IneqMeshViewModel ineqMeshViewModel = new IneqMeshViewModel 
            {
                IneqMesh = new MeshSamples()[Mesh],
                Quality = true,
                CurvatureQuality = true
            };
            Session["IneqMeshViewModel"] = ineqMeshViewModel;

            string path = System.IO.Path.Combine(Server.MapPath("~/Samples"), Mesh + ".ply");

            if (System.IO.File.Exists(path))
            {
                ineqMeshViewModel.PLY = System.IO.File.ReadAllText(path);
            }
            else
            {
                ineqMeshViewModel.CreateMesh();
                System.IO.File.WriteAllText(path, ineqMeshViewModel.PLY);
            }

            return Json(new { success = true }, JsonRequestBehavior.AllowGet);
        }

        public ActionResult SetSampleMesh1(string Mesh = "")
        {
            IneqMeshViewModel ineqMeshViewModel = new IneqMeshViewModel
            {
                //IneqMesh = new MeshSamples()[Mesh],
                //Quality = true,
                //CurvatureQuality = true
            };
            Session["IneqMeshViewModel"] = ineqMeshViewModel;

            string path = System.IO.Path.Combine(Server.MapPath("~/Samples"), Mesh + ".ply");

            if (System.IO.File.Exists(path))
            {
                ineqMeshViewModel.PLY = System.IO.File.ReadAllText(path);
                return Json(new { success = true }, JsonRequestBehavior.AllowGet);
            }
            else
            {
                ineqMeshViewModel.PLY = "";
                return Json(new { success = false }, JsonRequestBehavior.AllowGet);
            }            
        }

        [HttpPost, ValidateInput(false)]
        public ActionResult SetIneqMesh(IneqMeshViewModel ineqMeshViewModel, bool save = false)
        {
            if (ModelState.IsValid)
            {
                try
                {
                    if (save)
                    {
                        ((IneqMeshViewModel)Session["IneqMeshViewModel"]).Save(ineqMeshViewModel.Name);
                    }
                    else 
                    {
                        ineqMeshViewModel.SetIneqMesh();
                        Session["IneqMeshViewModel"] = ineqMeshViewModel;
                    }
                }
                catch (Exception exc)
                {
                    return Json(new { success = false, message = exc.Message });
                }

                return Json(new { success = true });
            }

            return Json(new { success = false, message = "Invalid data" });
        }

        public ActionResult GetMeshData(string name)
        {
            string path = System.IO.Path.Combine(Server.MapPath("~/Samples"), name + ".json");

            var json = System.IO.File.ReadAllText(path);

            return new ContentResult
            {
                Content = json,
                ContentType = "application/json",
                ContentEncoding = Encoding.UTF8
            };
        }

        public ActionResult GetMesh(bool stl = false)
        {
            IneqMeshViewModel ineqMeshViewModel = (IneqMeshViewModel)Session["IneqMeshViewModel"];

            if (ineqMeshViewModel == null)
            {
                return Content(null);
            }

            if (String.IsNullOrEmpty(ineqMeshViewModel.PLY))
            {
                ineqMeshViewModel.CreateMesh();
            }

            if (stl)
            {
                return File(PLYTools.ConvertPLYToSTL(ineqMeshViewModel.PLY), "application/octet-stream", "mesh.stl");
            }
            else
            {
                return File(System.Text.Encoding.UTF8.GetBytes(ineqMeshViewModel.PLY), "application/octet-stream", "mesh.ply");
            }
        }

        public ActionResult ImproveQuality(bool boundary)
        {
            IneqMeshViewModel ineqMeshViewModel = (IneqMeshViewModel)Session["IneqMeshViewModel"];

            if (ineqMeshViewModel == null || ineqMeshViewModel.IneqMesh == null)
            {
                return Content(null);
            }

            var ineqMesh = ineqMeshViewModel.IneqMesh;

            if (ineqMesh.Tetrahedrons.Count == 0) //bad check of mesh created
            {
                ineqMeshViewModel.CreateMesh(false);
            }

            if (boundary)
            {
                IneqMeshTools.CheckBoundaryQuality(ineqMesh);
            }
            else
            {
                IneqMeshTools.CheckQuality(ineqMesh);
            }

            ineqMeshViewModel.PLY = PLYTools.GetPLY(ineqMesh);

            return Content(null);
        }

        public ActionResult ImproveCurvatureQuality()
        {
            IneqMeshViewModel ineqMeshViewModel = (IneqMeshViewModel)Session["IneqMeshViewModel"];

            if (ineqMeshViewModel == null || ineqMeshViewModel.IneqMesh == null)
            {
                return Content(null);
            }

            var ineqMesh = ineqMeshViewModel.IneqMesh;

            if (ineqMesh.Tetrahedrons.Count == 0) //bad check of mesh created
            {
                ineqMeshViewModel.CreateMesh(false);
            }

            IneqMeshTools.CheckCurvatureQuality(ineqMesh);
            //IneqMeshTools.CheckCurvatureQuality(ineqMesh);            

            ineqMeshViewModel.PLY = PLYTools.GetPLY(ineqMesh);

            return Content(null);
        }


        public ActionResult Jiggle()
        {
            IneqMeshViewModel ineqMeshViewModel = (IneqMeshViewModel)Session["IneqMeshViewModel"];

            if (ineqMeshViewModel == null || ineqMeshViewModel.IneqMesh == null)
            {
                return Content(null);
            }

            var ineqMesh = ineqMeshViewModel.IneqMesh;

            if (ineqMesh.Tetrahedrons.Count == 0) //bad check of mesh created
            {
                ineqMeshViewModel.CreateMesh(false);
            }

            ineqMesh.Jiggle(3); 

            ineqMeshViewModel.PLY = PLYTools.GetPLY(ineqMesh);

            return Content(null);
        }

        public ActionResult Refine()
        {
            IneqMeshViewModel ineqMeshViewModel = (IneqMeshViewModel)Session["IneqMeshViewModel"];

            if (ineqMeshViewModel == null || ineqMeshViewModel.IneqMesh == null)
            {
                return Content(null);
            }

            var ineqMesh = ineqMeshViewModel.IneqMesh;

            if (ineqMesh.Tetrahedrons.Count == 0) //bad check of mesh created
            {
                ineqMeshViewModel.CreateMesh(false);
            }

            ineqMesh.RefineBoundaryTriangles(ineqMesh.Tetrahedrons.SelectMany(t => t.Triangles().Where(tr => tr.Boundary)));
            ineqMesh.DeleteLonelyPoints();
            ineqMesh.Jiggle(3, false);

            ineqMeshViewModel.PLY = PLYTools.GetPLY(ineqMesh);

            return Content(null);
        }

        //public ActionResult GetSampleFormula(int sampleIndex)
        //{
        //    var tmp = IneqMeshViewModel.DefaultModel(sampleIndex);

        //    return Json(new
        //    {
        //        formula = tmp.Formula
        //    });
        //}

        public ActionResult GetSampleUFunc(int sampleUFuncIndex)
        {
            var tmp = IneqMeshViewModel.DefaultModel(-1, sampleUFuncIndex);

            return Json(new
            {
                ufunc = tmp.UFunc
            });
        }
    }
}