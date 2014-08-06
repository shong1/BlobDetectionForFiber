#ifndef __itkBinaryTreeNodeBase_h
#define __itkBinaryTreeNodeBase_h

#include "itkTreeNode.h"

namespace itk
{
/** \class BinaryTreeNodeBase
 * \brief Base class for binary tree nodes, to be used with an herited class
 * from BinaryTreeContainer
 *
 * This class represents a binary tree node. All the mechanisms are in the node.
 * The Depth and Size values are calculated in O(1), because updated at very
 * shape change of the tree.
 *
 * Rotations are also implemented.
 *
 */
template< class TValue >
class BinaryTreeNodeBase : public TreeNode< TValue >
{
public:
	/** Standard typedefs. */
	typedef BinaryTreeNodeBase					Self;
	typedef SmartPointer< Self >				Pointer;
	typedef SmartPointer< const Self >			ConstPointer;
	typedef TreeNode< TValue >					Superclass;
	typedef Self								Baseclass;

	/** Method for creation through the object factory. */
	itkNewMacro( Self );

	/** Run-time type information (and related methods). */
	itkTypeMacro( BinaryTreeNodeBase, TreeNode );

	/** Convenient typedefs. */
	typedef TValue									ValueType;

	typedef typename Superclass::ChildIdentifier	ChildIdentifier;
	typedef typename Superclass::ChildrenListType	ChildrenListType;

	/** Right and Left indexes. */
	static const ChildIdentifier Left = 0, Right = 1;
	static const int MaxDepth = -1;

	/** Set/Get macros for the depth of the node, in O(1). */
	itkSetMacro( Depth, int );
	int GetDepth() const;

	/** Set/Get macros for the size of the subtree of the node, in O(1). */
	itkSetMacro( Size, int );
	int GetSize() const;

	/** Routine to add a child to the node. */
	void AddChild( ChildIdentifier position, Self *node );

	/** Right rotation. This must be called from a root node, and it
	 * returns the root, so it can be updated if it has been changed.
	 * The node parameter is the location where we do the rotation.
	 */
	template< class TSelf = Self >
	TSelf* RotateRight( TSelf *node );

	/** Right rotation. This must be called from a root node, and it
	 * returns the root, so it can be updated if it has been changed.
	 * The node parameter is the location where we do the rotation.
	 */
	template< class TSelf = Self >
	TSelf* RotateLeft( TSelf *node );

	/** Rotation. This must be called from a root node, and it
	 * returns the root, so it can be updated if it has been changed.
	 * The node parameter is the location where we do the rotation.
	 * The side parameter should be Right or Left
	 */
	template< class TSelf = Self >
	TSelf* Rotate( TSelf *node, ChildIdentifier side );

	/** Destroys all the children of the node (not the whole subtree). */
	virtual void WipeChildren();

	/** Returns the node storing elemenet. */
	template< class TSelf = Self >
	TSelf* GetNode( const ValueType element ) const;

	/** Returns the child of the node. The parameter should be Right or Left. */
	template< class TSelf = Self >
	TSelf * GetChild( ChildIdentifier number ) const { return dynamic_cast< TSelf* > ( Superclass::GetChild( number ) ); };

	/** Returns the parent of the node. If the node is a root, it returns a NULL pointer. */
	template< class TSelf = Self >
	TSelf * GetParent() const { return dynamic_cast< TSelf* > ( Superclass::GetParent() ); };

	/** Set the new parent. */
	void SetParent( Self * node );

	/** Routine to update the depth of all the ancesters node depth when
	 * the tree is changed.
	 */
	virtual void UpdateDepthToRoot( Self *node );

	/** Increments the depth of the node. */
	virtual void IncrementDepth();

	/** Decrements the depth of the node. */
	virtual void DecrementDepth();

	/** Add a node to the tree. It should be called from the root.
	 * It returns the new root. Not implemented because it is a
	 * base class. It returns NULL.
	 */
	template< class TSelf = Self >
	TSelf* Add( Self *node ){ return NULL; }

	/** Remove a node to the tree. It should be called from the root.
	 * It returns the new root. Not implemented because it is a
	 * base class. It returns NULL.
	 */
	template< class TSelf = Self >
	TSelf* Remove( Self *node ){ return NULL; }

	/** Splits the tree at node location. It should be called from the root.
	 * It returns the new root. Not implemented because it is a
	 * base class. It returns NULL.
	 * The pointer node is in the right side, and root gets the left one.
	 */
	template< class TSelf = Self >
	TSelf* Split( TSelf* node, TSelf **root ){ return NULL; }

	/** Merges with root a node to the tree. It should be called from the root.
	 * It returns the new root. Not implemented because it is a
	 * base class. It returns NULL.
	 */
	template< class TSelf = Self >
	TSelf* Merge( TSelf *root ){ return NULL; }

	/** Shares node's information to self (depth, size, etc.). */
	virtual void ShareInfo( Self * node ){}

	/** Copies other element. */
	virtual void Set( const Self & other );

	/** Set element to self. */
	virtual void Set( const ValueType & element );

	template< class TSelf = Self >
	bool operator<( const TSelf &other );
	template< class TSelf = Self >
	bool operator>( const TSelf &other );
	template< class TSelf = Self >
	bool operator<=( const TSelf &other );
	template< class TSelf = Self >
	bool operator>=( const TSelf &other );
	template< class TSelf = Self >
	bool operator==( const TSelf &other );
	template< class TSelf = Self >
	bool operator!=( const TSelf &other );
	template< class TSelf = Self >
	bool operator<( const ValueType & element );
	template< class TSelf = Self >
	bool operator>( const ValueType & element );
	template< class TSelf = Self >
	bool operator<=( const ValueType & element );
	template< class TSelf = Self >
	bool operator>=( const ValueType & element );
	template< class TSelf = Self >
	bool operator==( const ValueType & element );
	template< class TSelf = Self >
	bool operator!=( const ValueType & element );

protected:
	BinaryTreeNodeBase();
	~BinaryTreeNodeBase(){};
	int m_Depth, m_Size;

private:
	BinaryTreeNodeBase( const Self & ); // purposely node implemented
	void operator=( const Self & ); // purposely node implemented

};

}// namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBinaryTreeNodeBase.hxx"
#endif

#endif
