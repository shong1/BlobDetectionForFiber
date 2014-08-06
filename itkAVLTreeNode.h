#ifndef __itkAVLTreeNode_h
#define __itkAVLTreeNode_h

#include "itkSearchBinaryTreeNode.h"

namespace itk
{
/** \class AVLTreeNode
 * \brief AVL tree node implementation, to be used with an herited class
 * from BinaryTreeContainer
 *
 * This class represents an AVL tree node. All the mechanisms are in the node.
 * The Depth and Size values are calculated in O(1), because updated at very
 * shape change of the tree. An AVL tree is an auto balanced tree.
 *
 */
template< class TValue >
class AVLTreeNode : public SearchBinaryTreeNode< TValue >
{
public:
	/** Standard typedefs. */
	typedef AVLTreeNode								Self;
	typedef SmartPointer< Self >					Pointer;
	typedef SmartPointer< const Self >				ConstPointer;
	typedef SearchBinaryTreeNode< TValue >			Superclass;

	/** Method for creation through the object factory. */
	itkNewMacro( Self );

	/** Run-time type information (and related methods). */
	itkTypeMacro( AVLTreeNode, SearchBinaryTreeNode );

	/** Convenient typedefs. */
	typedef TValue									ValueType;

	typedef typename Superclass::ChildIdentifier	ChildIdentifier;
	typedef typename Superclass::ChildrenListType	ChildrenListType;
	typedef typename Superclass::Baseclass			Baseclass;

	/** Right and Left indexes. */
	static const ChildIdentifier Left = Superclass::Left, Right = Superclass::Right;
	static const int MaxDepth = Superclass::MaxDepth;

	/** Right rotation. This must be called from a root node, and it
	 * returns the root, so it can be updated if it has been changed.
	 * The node parameter is the location where we do the rotation.
	 */
	template< class TSelf = Self >
	TSelf* RotateRight( TSelf *node ){ return Superclass::template RotateRight< TSelf >( node ); }

	/** Right rotation. This must be called from a root node, and it
	 * returns the root, so it can be updated if it has been changed.
	 * The node parameter is the location where we do the rotation.
	 */
	template< class TSelf = Self >
	TSelf* RotateLeft( TSelf *node ){ return Superclass::template RotateLeft< TSelf >( node ); }

	/** Equilibrate the tree from a node. Should be called from the root.
	 * It returns the new root in case it changed. The recursive parameter
	 * tells if we equilibrate the whole tree, or just the first unbalanced
	 * node.
	 */
	template< class TSelf = Self >
	TSelf* Equilibrate( TSelf *node, bool recursive = false );

	/** Returns the node storing elemenet. */
	template< class TSelf = Self >
	TSelf* GetNode( const ValueType element ){ return Superclass::template GetNode< TSelf >( element ); }

	/** Returns a pointer to the node lower of even to element. */
	template< class TSelf = Self >
	TSelf* GetLowerNode( const ValueType element ){ return Superclass::template GetLowerNode< TSelf >( element ); }

	/** Returns a pointer to the node upper of even to element. */
	template< class TSelf = Self >
	TSelf* GetUpperNode( const ValueType element ){ return Superclass::template GetUpperNode< TSelf >( element ); }

	/** Returns a pointer to the leaf strictly lower to element. */
	template< class TSelf = Self >
	TSelf* GetStrictlyLowerNode( const ValueType element ){ return Superclass::template GetStrictlyLowerNode< TSelf >( element ); }

	/** Returns a pointer to the leaf strictly upper to element. */
	template< class TSelf = Self >
	TSelf* GetStrictlyUpperNode( const ValueType element ){ return Superclass::template GetStrictlyUpperNode< TSelf >( element ); }

	/** Return the lowest node of a certain depth. */
	template< class TSelf = Self >
	TSelf* GetLowestNode( int depth = MaxDepth ){ return Superclass::template GetLowestNode< TSelf >( depth ); }

	/** Return the uppest node of a certain depth. */
	template< class TSelf = Self >
	TSelf* GetUppestNode( int depth = MaxDepth ){ return Superclass::template GetUppestNode< TSelf >( depth ); }

	/** Returns a pointer to the leaf lower of even to element. */
	template< class TSelf = Self >
	TSelf* GetLowerLeaf( const ValueType element ){ return Superclass::template GetLowerLeaf< TSelf >( element ); }

	/** Returns a pointer to the leaf upper or even to element. */
	template< class TSelf = Self >
	TSelf* GetUpperLeaf( const ValueType element ){ return Superclass::template GetUpperLeaf< TSelf >( element ); }

	/** Returns a pointer to the leaf strictly lower to element. */
	template< class TSelf = Self >
	TSelf* GetStrictlyLowerLeaf( const ValueType element ){ return Superclass::template GetStrictlyLowerLeaf< TSelf >( element ); }

	/** Returns a pointer to the leaf strictly upper to element. */
	template< class TSelf = Self >
	TSelf* GetStrictlyUpperLeaf( const ValueType element ){ return Superclass::template GetStrictlyUpperLeaf< TSelf >( element ); }

	/** Returns the parent of the node. If the node is a root, it returns a NULL pointer. */
	template< class TSelf = Self >
	TSelf* GetParent() const { return Superclass::template GetParent< TSelf >(); }

	/** Returns the child of the node. The parameter should be Right or Left. */
	template< class TSelf = Self >
	TSelf* GetChild( ChildIdentifier number ) const { return Superclass::template GetChild< TSelf >( number ); }

	/** Shares node's information to itself. */
	virtual void ShareInfo( Self* node ){}

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

	/** Merges with root a node to the tree. It should be called from the root.
	 * It returns the new root.
	 */
	template< class TSelf = Self >
	TSelf* Merge( TSelf *root );

	/** Splits the tree at node location. It should be called from the root.
	 * It returns the new root.
	 * The pointer node is in the right side, and root gets the left one.
	 */
	template< class TSelf = Self >
	TSelf* Split( TSelf *node, TSelf **root ){ return Superclass::template Split< TSelf >( node, root ); }


protected:
	AVLTreeNode();
	~AVLTreeNode(){};

private:
	AVLTreeNode( const Self & ); // purposely node implemented
	void operator=( const Self & ); // purposely node implemented

};

}// namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAVLTreeNode.hxx"
#endif

#endif
