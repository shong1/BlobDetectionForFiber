#ifndef __itkHalfConvexHullTreeNode_h
#define __itkHalfConvexHullTreeNode_h

#include "itkAVLVectorIntervalTreeNode.h"	
#include "itkHalfConvexHullTreeContainer.h"

namespace itk
{
/** \class HalfConvexHullTreeNode
 * \brief Node structure for computation of the half of a convex hull in \f$ O(\log^2n) \f$.
 *
 * This class uses an AVL node structure to compute the half of a convex hull
 * as described in the paper Maintenance of Configurations in the Plane, by
 * Mark H. OverMars and Jan Van Leeuwen :
 * http://www.ics.uci.edu/~goodrich/teach/cs164/notes/sdarticle-37.pdf
 *
 * It computes dynamically the insertion and suppression of a point
 * in \f$ O(\log^2n)\f$ in the worst case.
 *
 */
template< class TVector >
class HalfConvexHullTreeNode : public AVLTreeNode< TVector >
{
public:
	/** Standard typedefs. */
	typedef HalfConvexHullTreeNode							Self;
	typedef SmartPointer< Self >							Pointer;
	typedef SmartPointer< const Self >						ConstPointer;
	typedef AVLTreeNode< TVector >							Superclass;

	/** Method for creation through the object factory. */
	itkNewMacro( Self );

	/** Run-time type information (and related methods). */
	itkTypeMacro( HalfConvexHullTreeNode, AVLTreeNode );

	/** Convenient typedefs. */
	typedef TVector											VectorType;
	typedef TVector											ValueType;
	typedef typename TVector::ValueType						ScalarType;

	typedef typename Superclass::ChildIdentifier			ChildIdentifier;
	typedef typename Superclass::ChildrenListType			ChildrenListType;
	typedef typename Superclass::Baseclass					Baseclass;

	typedef BinaryTreeContainer< AVLVectorIntervalTreeNode< TVector > > ConcatenableQueueType;
	typedef typename ConcatenableQueueType::TreeNodeType	ConcatenableQueueNodeType;

	/** Right and Left indexes. */
	static const ChildIdentifier Left = Superclass::Left, Right = Superclass::Right;

	/** Get the convex hull. */
	ConcatenableQueueType* GetQ() const { return m_Q; }

	/** Set the convex hull. */
	void SetQ( ConcatenableQueueType* q ) { m_Q = q; }

	/** Get the splitting bridge for maintenance. See paper. */
	ConcatenableQueueNodeType* GetBridge() const { return m_Bridge; }

	/** Set the splitting bridge. See paper. */
	void SetBridge( ConcatenableQueueNodeType* bridge ) { m_Bridge = bridge; }

	/** Get the whole tree to get some global information. */
	const HalfConvexHullTreeContainer< Self >* GetTree() const { return m_Tree; }

	/** Set the tree to have access to some global information. */
	void SetTree( const HalfConvexHullTreeContainer< Self >* tree ) { m_Tree = tree; }
	template< class TConcatenableQueueNode = ConcatenableQueueNodeType >
	TConcatenableQueueNode* GetQRoot() const { return dynamic_cast< TConcatenableQueueNode* >( this->GetQ()->GetRoot() ); }

	/** Right rotation. This must be called from a root node, and it
	 * returns the root, so it can be updated if it has been changed.
	 * The node parameter is the location where we do the rotation.
	 */
	template< class TSelf = Self >
	TSelf* RotateRight( TSelf *node );
	/** Left rotation. This must be called from a root node, and it
	 * returns the root, so it can be updated if it has been changed.
	 * The node parameter is the location where we do the rotation.
	 */
	template< class TSelf = Self >
	TSelf* RotateLeft( TSelf *node );

	/** Equilibrate the tree from a node. Should be called from the root.
	 * It returns the new root in case it changed.
	 */
	template< class TSelf = Self >
	TSelf* Equilibrate( TSelf *node );

	/** Returns the child of the node. The parameter should be Right or Left. */
	template< class TSelf = Self >
	TSelf * GetChild( ChildIdentifier number ) const { return Superclass::template GetChild< TSelf >( number ); }

	/** Returns the parent of the node. If the node is root, it returns a NULL pointer. */
	template< class TSelf = Self >
	TSelf * GetParent() const { return Superclass::template GetParent< TSelf >(); };

	/** Add a node to the hull. It should be called from the root.
	 * It returns the new root.
	 */
	template< class TSelf = Self >
	TSelf * Add( ConcatenableQueueNodeType * node, TSelf* newRoot = NULL );

	/** Remove a node to the hull. It should be called from the root.
	 * It returns the new root.
	 */
	template< class TSelf = Self >
	TSelf * Remove( const VectorType & element );

	/** Down routine. See paper. */
	template< class TSelf = Self >
	TSelf * Down( TSelf * to );

	/** Up routine. See paper. */
	template< class TSelf = Self >
	TSelf * Up( TSelf * node, int number = -1 );

	/** Tells if the point x is "under" the vector (down, up) or "above". It returns
	 * Right for "above" and Left for "under".
	 */
	bool ContourSide( const VectorType & down, const VectorType & up, const VectorType & x, ScalarType * result ) const;

	/** Tells if the chidren of the current node have no queue inside. */
	bool ChildrenAreEmpty() const;

	/** Routine to update the bridge. See paper. */
	void UpdateBridge( ConcatenableQueueNodeType *i, ConcatenableQueueNodeType *j, ConcatenableQueueNodeType **bridge, ConcatenableQueueNodeType **brjdge );

protected:
	HalfConvexHullTreeNode(){ m_Bridge = NULL; }
	~HalfConvexHullTreeNode(){};

	/** Hull of the current node. */
	SmartPointer< ConcatenableQueueType > m_Q;

	/** Pointer to the bridge. */
	ConcatenableQueueNodeType *m_Bridge;

	/** Pointer to the main tree container. */
	const HalfConvexHullTreeContainer< Self > *m_Tree;
	

private:
	HalfConvexHullTreeNode( const Self & ); // purposely node implemented
	void operator=( const Self & ); // purposely node implemented

};

}// namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkHalfConvexHullTreeNode.hxx"
#endif

#endif
