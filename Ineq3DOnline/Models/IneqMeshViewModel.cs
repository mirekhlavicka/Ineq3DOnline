using MeshData;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Web;

namespace Ineq3DOnline.Models
{
    public class IneqMeshViewModel
    {
        public string Formula { get; set; }

        public bool Quality { get; set; }

        public IneqMesh IneqMesh  = null;

        public string PLY { get; set; }
    }
}