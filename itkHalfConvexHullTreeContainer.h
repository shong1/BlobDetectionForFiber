#ifndef __itkHalfConvexHullTreeContainer_h
#define __itkHalfConvexHullTreeContainer_h

#include "itkTreeContainerBase.h"
#include "itkBinaryTreeContainer.h"
#include "itkAVLVectorIntervalTreeNode.h"
#include "itkTreeNode.h"

namespace itk
{
template< class TTreeNode >
class HalfConvexHullTreeContainer : public TreeContainerBase< typename TTreeNode::ValueType >
{
/** \class HalfConvexHullTreeContainer
 * \brief Main strucure to use with an herited class of HalfConvexHullTreeNode
 * to have a dynamically half convex hull updated in \f$O(\log^2n)\f$.
 *
 * This class is the main structure using HalfConvexHullTreeNode nodes.
 * It can either add or remove a point in 2D.
 * m_Dimension and m_OtherDimension allow us to use n dimensional vector,
 * saying in which projection we want to work.
 *
 */
public:
	/** Standard typedefs */
	typedef HalfConvexHullTreeContainer									Self;
	typedef SmartPointer< Self >								Pointer;
	typedef SmartPointer< const Self >							ConstPointer;
	typedef TreeContainerBase< typename TTreeNode::ValueType >	Superclass;

	/** Method for creation through the object factory. */
	itkNewMacro( Self );

	/** Run-time type information (and related methods). */
	itkTypeMacro( HalfConvexHullTreeContainer, TreeContainerBase );

	/** Convenient typedefs. */
	typedef TTreeNode								TreeNodeType;
	typedef typename TreeNodeType::ValueType		ValueType;
	typedef ValueType								VectorType;
	typedef typename TreeNodeType::Pointer			TreeNodePointer;
	typedef typename TreeNodeType::ChildIdentifier	ChildIdentifier;
	typedef TreeNode< ValueType >					TreeNodeBaseType;
	typedef AVLVectorIntervalTreeNode< VectorType >	ConcatenableQueueNodeType;
	typedef BinaryTreeContainer< ConcatenableQueueNodeType >		ConcatenableQueueType;

	/** Right and Left Indexes. */
	static const ChildIdentifier Right = TreeNodeType::Right, Left = TreeNodeType::Left;

	/** Dimension of the vector used. */
	static const unsigned int Dimension = VectorType::Dimension;

	/** Gets rid of every nodes of the tree. */
	virtual bool Clear();
	virtual bool Clear( TreeNodeBaseType * node );

	/** Tells if the tree contains an element. */
	virtual bool Contains( const ValueType element );

	/** Counts the number of nodes in the tree. */
	virtual int Count() const;
	virtual int Count( TreeNodeType * node ) const;

	/** Returns the root of the trees. */
	TreeNodeBaseType * GetRoot() const;
	
	/** Tells if the element is a leaf. */
	virtual bool IsLeaf( const ValueType element );

	/** Tells if the element is a root. */
	virtual bool IsRoot( const ValueType element );

	/** Sets the root of the tree. */
	virtual bool SetRoot( const ValueType element );
	virtual bool SetRoot( TreeNodeBaseType * node );

	/** Set/Get macros to know the side of the contour.
	 * It should be Right or Left.
	 */
    itkGetConstMacro( Contour, ChildIdentifier );
	itkSetMacro( Contour, ChildIdentifier );

	/** Set/Get macros to know in which dimension we work
	 * to sort the nodes.
	 */
	itkGetConstMacro( Dimension, unsigned int );
	itkSetMacro( Dimension, unsigned int );

	/** Set/Get macros to know which is the "depth" dimension. */
	itkGetConstMacro( OtherDimension, unsigned int );
	itkSetMacro( OtherDimension, unsigned int );

	/** Get the root of the convex hull (root node). */
	ConcatenableQueueNodeType* GetQRoot();

	/** Add an element to the set of points. */
	virtual ConcatenableQueueNodeType* Add( const VectorType element );

	/** Add a node to the set of points. */
	virtual void Add( ConcatenableQueueNodeType *node );

	/** Remove an element of the set of points. */
	virtual ConcatenableQueueNodeType* Remove( const VectorType element );


protected:
	HalfConvexHullTreeContainer(){ m_Root = NULL; };
	virtual ~HalfConvexHullTreeContainer(){};

	/** Root of the tree. */
	TreeNodePointer m_Root;

	/** Contour side for the half convex hull. */
	ChildIdentifier m_Contour;

	/** Dimensions to work with. */
	unsigned int m_Dimension, m_OtherDimension;

private:
	HalfConvexHullTreeContainer( const Self & );	// purposely not implemented
	void operator=( const Self & );			// purposely not implemented

};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkHalfConvexHullTreeContainer.hxx"
#endif

#endif
