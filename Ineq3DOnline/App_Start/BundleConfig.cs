using System.Web;
using System.Web.Optimization;

namespace Ineq3DOnline
{
    public class BundleConfig
    {
        // For more information on bundling, visit http://go.microsoft.com/fwlink/?LinkId=301862
        public static void RegisterBundles(BundleCollection bundles)
        {
            bundles.Add(new ScriptBundle("~/bundles/jquery").Include(
                        "~/Scripts/jquery-{version}.js"));

            bundles.Add(new ScriptBundle("~/bundles/jqueryval").Include(
                        "~/Scripts/jquery.validate*"));

            // Use the development version of Modernizr to develop with and learn from. Then, when you're
            // ready for production, use the build tool at http://modernizr.com to pick only the tests you need.
            bundles.Add(new ScriptBundle("~/bundles/modernizr").Include(
                        "~/Scripts/modernizr-*"));

            bundles.Add(new ScriptBundle("~/bundles/bootstrap").Include(
                      "~/Scripts/bootstrap.js",
                      "~/Scripts/respond.js"));

            bundles.Add(new ScriptBundle("~/bundles/parser").Include(
                      "~/Scripts/parser.js"));


            bundles.Add(new StyleBundle("~/Content/css").Include(
                      "~/Content/bootstrap.css",
                      "~/Scripts/codemirror/lib/codemirror.css",
                      "~/Scripts/codemirror/addon/fold/foldgutter.css",
                      "~/Content/site.css"));

            bundles.Add(new ScriptBundle("~/bundles/three").Include(
                        "~/Scripts/three/three.js",
                        "~/Scripts/three/TrackballControls.js",
                        "~/Scripts/three/CanvasRenderer.js",
                        "~/Scripts/three/Projector.js",
                        "~/Scripts/Spin.js",
                        "~/Scripts/three/PLYLoader.js"));

            bundles.Add(new ScriptBundle("~/bundles/codemirror")
                .Include("~/Scripts/codemirror/lib/codemirror.js")
                .Include("~/Scripts/codemirror/mode/clike/clike.js")
                .Include("~/Scripts/codemirror/addon/edit/matchbrackets.js")
                .Include("~/Scripts/codemirror/addon/edit/closebrackets.js")
                .Include("~/Scripts/codemirror/addon/fold/foldcode.js")
                .Include("~/Scripts/codemirror/addon/fold/foldgutter.js")
                .Include("~/Scripts/codemirror/addon/fold/brace-fold.js")
                .Include("~/Scripts/codemirror/addon/fold/indent-fold.js")
                .Include("~/Scripts/codemirror/addon/fold/comment-fold.js")
                );

        }
    }
}
