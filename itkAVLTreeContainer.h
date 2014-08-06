#ifndef __itkAVLTreeContainer_h
#define __itkAVLTreeContainer_h

#include "itkBinaryTreeContainer.h"

namespace itk
{
template< class TTreeNode >
class AVLTreeContainer : public BinaryTreeContainer< TTreeNode >
{
public:
	/* Standard typedefs */
	typedef AVLTreeContainer							Self;
	typedef SmartPointer< Self >						Pointer;
	typedef SmartPointer< const Self >					ConstPointer;
	typedef BinaryTreeContainer< TTreeNode >			Superclass;

	/* Mehtod for creation through the object factory. */
	itkNewMacro( Self );

	/* Run-time type information (and related methods). */
	itkTypeMacro( AVLTreeContainer, BinaryTreeContainer );

	typedef typename Superclass::TreeNodeType			TreeNodeType;
	typedef typename Superclass::TreeNodeBaseType		TreeNodeBaseType;
	typedef typename Superclass::TreeNodePointer		TreeNodePointer;
	typedef typename Superclass::ChildIdentifier		ChildIdentifier;
	
	typedef typename Superclass::ValueType				ValueType;

	static const ChildIdentifier Right = TreeNodeType::Right, Left = TreeNodeType::Left;

	virtual TreeNodeType * Add( const ValueType element );
	virtual void Add( TreeNodeType * node );

	virtual bool Equilibrate( TreeNodeType * node, bool recursive = false );
	virtual bool Equilibrate( TreeNodeType * node, TreeNodeType * root, bool recursive = false );

	virtual bool Remove( TreeNodeType * X );
	virtual bool Remove( const ValueType element );

	virtual void Merge( Self * tree );
	virtual void Merge( TreeNodeType * tree );
	Self* Split( const ValueType element , Self* tree );
	Self* SizeSplit( int size , Self* tree );

	TreeNodeType * GetRoot() { return dynamic_cast< TreeNodeType * >( Superclass::GetRoot() ); }

protected:
	AVLTreeContainer(){};
	virtual ~AVLTreeContainer(){};

private:
	AVLTreeContainer( const Self & ); // purposely not implemented
	void operator=( const Self & ); // purposely not implemented
};
}// namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAVLTreeContainer.hxx"
#endif

#endif
