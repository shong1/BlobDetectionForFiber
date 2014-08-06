#ifndef __itkConvexHullAugmentedTreeNode_h
#define __itkConvexHullAugmentedTreeNode_h

#include "itkAugmentedTreeNode.h"

namespace itk
{
template< class TBinaryTreeNode >
class ConvexHullAugmentedTreeNode : AugmentedTreeNode< TBinaryTreeNode >
{
public:
	typedef ConvexHullAugmentedTreeNode					Self;
	typedef SmartPointer< Self >						Pointer;
	typedef SmartPointer< const Self >					ConstPointer;
	typedef TBinaryTreeNode								Superclass;

	itkNewMacro( Self );

	itkTypeMacro( ConvexHullAugmentedTreeNode, AugmentedTreeNode );

	typedef typename Superclass::ValueType				ValueType;
	typedef typename Superclass::ChildIdentifier		ChildIdentifier;
	typedef typename Superclass::NodeContainerType		NodeContainerType;

	static const ChildIdentifier Left = Superclass::Left, Right = Superclass::Right;

	Self* GetSubtree( ChildIdentifier number ) const { return dynamic_cast< Self* >( Superclass::GetSubtree( number ) ); };
	Self* GetChild( ChildIdentifier number ) const { return dynamic_cast< Self* >( Superclass::GetChild( number ) ); };
	Self* GetParent() const { return dynamic_cast< Self* >( Superclass::GetParent() ); };
	Self* GetSubtreeParent() const { return dynamic_cast< Self* >( Superclass::GetSubtreeParent() ); };

	virutal void AddSubtree( ChildIdentifier number, Self * node );

protected:
	ConvexHullAugmentedTreeNode(){};
	~ConvexHullAugmentedTreeNode(){};

};
}// namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itlConvexHullAugmentedTreeNode.hxx"
#endif

#endif
