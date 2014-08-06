#ifndef __itkConvexHullAugmentedTreeNode_hxx
#define __itkConvexHullAugmentedTreeNode_hxx

#include "itkConvexHullAugmentedTreeNode.h"

namespace itk
{
template< class TBinaryTreeNode >
void ConvexHullAugmentedTreeNode< TBinaryTreeNode >
::AddSubtree( ChildIdentifier number, Self * node )
{
	if( number > 1 )
	{
		std::cout << "Error: itkConvexHullAugmentedTreeNode, AddSubtree( number, node ) : number > 1 );
	}
	Superclass::AddSubtree( number, node );
}

}// namespace itk

#endif
