using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace MeshData
{
    using FuncXYZ = Func<double, double, double, double>;

    public class IneqTree
    {
        private IneqNode root = null;
        private List<FuncXYZ> expressionList = null;

        public List<FuncXYZ> ExpressionList
        {
            get 
            {
                if (expressionList == null)
                {
                    expressionList = new List<FuncXYZ>();
                    AddToexpressionList(root);
                }

                return expressionList; 
            }
        }

        public IneqNode Root
        {
            get { return root; }
        }

        public IneqTree()
        { 
        }

        public IneqTree(NodeType nodeType, IneqTree left, IneqTree right)
        {
            if (left.root == null)
            {
                root = right.root;
            }
            else if (right.root == null)
            {
                root = left.root;
            }
            else
            {
                root = new IneqNode();
                root.Expression = null;
                root.NodeType = nodeType;
                root.Left = left.root;
                root.Right = right.root;
            }
        }

        public IneqTree(FuncXYZ expression)
        {
            root = new IneqNode();
            root.Expression = expression;
            root.NodeType = NodeType.NodeExpression;
            root.Left = null;
            root.Right = null;
        }

        private void AddToexpressionList(IneqNode node)
        {
            if (node == null)
                return;

            if (node.NodeType == NodeType.NodeExpression)
            {
                expressionList.Add(node.Expression);
                node.ExpressionIndex = expressionList.Count - 1;
            }

            AddToexpressionList(node.Left);
            AddToexpressionList(node.Right);

        }

        public static IneqTree operator &(IneqTree left, IneqTree right)
        {
            return new IneqTree(NodeType.NodeAnd, left, right);
        }

        public static IneqTree operator |(IneqTree left, IneqTree right)
        {
            return new IneqTree(NodeType.NodeOr, left, right);
        }

        public static IneqTree operator &(IneqTree left, FuncXYZ right)
        {
            return new IneqTree(NodeType.NodeAnd, left, right);
        }

        public static IneqTree operator &(FuncXYZ left, IneqTree right)
        {
            return new IneqTree(NodeType.NodeAnd, left, right);
        }

        public static IneqTree operator |(IneqTree left, FuncXYZ right)
        {
            return new IneqTree(NodeType.NodeOr, left, right);
        }

        public static IneqTree operator |(FuncXYZ left, IneqTree right)
        {
            return new IneqTree(NodeType.NodeOr, left, right);
        }

        public static implicit operator IneqTree(FuncXYZ expression)
        {
            return new IneqTree(expression);
        }
        
        public enum NodeType
        {
            NodeExpression,
            NodeOr,
            NodeAnd
        }

        public class IneqNode
        {
            private NodeType nodeType;
            private FuncXYZ expression = null;
            private int expressionIndex = -1;
            private IneqNode left = null;
            private IneqNode right = null;

            public NodeType NodeType
            {
                get { return nodeType; }
                set { nodeType = value; }
            }

            public FuncXYZ Expression
            {
                get { return expression; }
                set { expression = value; }
            }

            public int ExpressionIndex
            {
                get { return expressionIndex; }
                set { expressionIndex = value; }
            }

            public IneqNode Left
            {
                get { return left; }
                set { left = value; }
            }

            public IneqNode Right
            {
                get { return right; }
                set { right = value; }
            }

            public IEnumerable<int> ExpressionIndexes
            {
                get
                {
                    if (nodeType == IneqTree.NodeType.NodeExpression)
                        return new int[] { expressionIndex };
                    else
                        return left.ExpressionIndexes.Union(right.ExpressionIndexes);
                }
            }

            public int SubNodeCount
            {
                get
                {
                    if (nodeType == IneqTree.NodeType.NodeExpression)
                        return 0;
                    else
                        return left.SubNodeCount + right.SubNodeCount;
                
                }
            }
        }
    }
}
