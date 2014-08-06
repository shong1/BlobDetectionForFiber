#ifndef __itkSearchBinaryTreeNode_h
#define __itkSearchBinaryTreeNode_h

#include "itkBinaryTreeNodeBase.h"

namespace itk
{
/** \class SearchBinaryTreeNode
 * \brief Binary search tree node implementation. To be used with an herited
 * class from BinaryTreeContainer
 *
 * This class has binary search trees implementations, and have thus algorithms
 * in \f$ O(\log(n)) \f$.
 *
 */
template< class TValue >
class SearchBinaryTreeNode : public BinaryTreeNodeBase< TValue >
{
public:
	/** Standard typedefs. */
	typedef SearchBinaryTreeNode					Self;
	typedef SmartPointer< Self >					Pointer;
	typedef SmartPointer< const Self >				ConstPointer;
	typedef BinaryTreeNodeBase< TValue >			Superclass;

	/** Method for creation through the object factory. */
	itkNewMacro( Self );

	/** Run-time type information (and related methods). */
	itkTypeMacro( SearchBinaryTreeNode, BinaryTreeNodeBase );

	/** Convenient typedefs. */
	typedef TValue									ValueType;

	typedef typename Superclass::ChildIdentifier	ChildIdentifier;
	typedef typename Superclass::ChildrenListType	ChildrenListType;
	typedef typename Superclass::Baseclass			Baseclass;

	static const ChildIdentifier Left = Superclass::Left, Right = Superclass::Right;
	static const int MaxDepth = Superclass::MaxDepth;

	/** Returns a pointer to the leaf lower of even to element. */
	template< class TSelf = Self >
	TSelf* GetLowerLeaf( const ValueType element );

	/** Returns a pointer to the leaf upper or even to element. */
	template< class TSelf = Self >
	TSelf* GetUpperLeaf( const ValueType element );

	/** Returns a pointer to the leaf strictly lower to element. */
	template< class TSelf = Self >
	TSelf* GetStrictlyLowerLeaf( const ValueType element );

	/** Returns a pointer to the leaf strictly upper to element. */
	template< class TSelf = Self >
	TSelf* GetStrictlyUpperLeaf( const ValueType element );

	/** Returns a pointer to the node lower of even to element. */
	template< class TSelf = Self >
	TSelf* GetLowerNode( const ValueType element );

	/** Returns a pointer to the node upper of even to element. */
	template< class TSelf = Self >
	TSelf* GetUpperNode( const ValueType element );

	/** Returns a pointer to the leaf strictly lower to element. */
	template< class TSelf = Self >
	TSelf* GetStrictlyLowerNode( const ValueType element );

	/** Returns a pointer to the leaf strictly upper to element. */
	template< class TSelf = Self >
	TSelf* GetStrictlyUpperNode( const ValueType element );

	/** Return the lowest node of a certain depth. */
	template< class TSelf = Self >
	TSelf* GetLowestNode( int depth = MaxDepth );

	/** Return the uppest node of a certain depth. */
	template< class TSelf = Self >
	TSelf* GetUppestNode( int depth = MaxDepth );

	/** Returns the node storing elemenet. */
	template< class TSelf = Self >
	TSelf* GetNode( const ValueType element );

	/** Shares node's information to self (depth, size, etc.). */
	virtual void ShareInfo( Self * node ){}

	/** Returns the parent of the node. If the node is a root, it returns a NULL pointer. */
	template< class TSelf = Self >
	TSelf* GetParent() const { return Superclass::template GetParent< TSelf >(); }

	/** Returns the child of the node. The parameter should be Right or Left. */
	template< class TSelf = Self >
	TSelf* GetChild( ChildIdentifier number ) const { return Superclass::template GetChild< TSelf >( number ); }

	/** Add a node to the tree. It should be called from the root.
	 * It returns the new root.
	 */
	template< class TSelf = Self >
	TSelf* Add( TSelf *node );

	/** Remove a node to the tree. It should be called from the root.
	 * It returns the new root.
	 */
	template< class TSelf = Self >
	TSelf* Remove( TSelf *node );

	/** Splits the tree at node location. It should be called from the root.
	 * It returns the new root.
	 * The pointer node is in the right side, and root gets the left one.
	 */
	template< class TSelf = Self >
	TSelf* Split( TSelf *node, TSelf **root );

	/** Merges with root a node to the tree. It should be called from the root.
	 * It returns the new root.
	 */
	template< class TSelf = Self >
	TSelf* Merge( TSelf *root );

protected:
	SearchBinaryTreeNode();
	~SearchBinaryTreeNode(){};

private:
	SearchBinaryTreeNode( const Self & ); // purposely node implemented
	void operator=( const Self & ); // purposely node implemented

};

}// namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSearchBinaryTreeNode.hxx"
#endif

#endif
