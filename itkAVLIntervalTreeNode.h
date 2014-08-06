#ifndef __itkAVLIntervalTreeNode_h
#define __itkAVLIntervalTreeNode_h

#include "itkAVLTreeNode.h"

namespace itk
{
/** \class AVLIntervalTreeNode
 * \brief AVL interval tree node implementation, to be used with an herited class
 * from BinaryTreeContainer
 *
 * This class represents an AVL interval tree node. All the mechanisms are in the node.
 * The Depth and Size values are calculated in O(1), because updated at very
 * shape change of the tree. An AVL tree is an auto balanced tree.
 *
 * The structure of this tree is as follow : every "elements" are at the leaves
 * (they are singletons). And every interval node's subtree is composed of nodes
 * that are in the root interval.
 *
 */
template< class TValue >
class AVLIntervalTreeNode : public AVLTreeNode< TValue >
{
public:
	/** Standard typedefs. */
	typedef AVLIntervalTreeNode						Self;
	typedef SmartPointer< Self >					Pointer;
	typedef SmartPointer< const Self >				ConstPointer;
	typedef AVLTreeNode< TValue >					Superclass;

	/** Method for creation through the object factory. */
	itkNewMacro( Self );

	/** Run-time type information (and related methods). */
	itkTypeMacro( AVLIntervalTreeNode, AVLTreeNode );

	/** Convenient typedefs. */
	typedef TValue									ValueType;

	typedef typename Superclass::ChildIdentifier	ChildIdentifier;
	typedef typename Superclass::ChildrenListType	ChildrenListType;
	typedef typename Superclass::Baseclass			Baseclass;

	/** Right and Left indexes. */
	static const ChildIdentifier Left = Superclass::Left, Right = Superclass::Right;
	static const int MaxDepth = Superclass::MaxDepth;

	/** Get the lower value of the interval. */
	virtual ValueType GetLow() const { return m_Low->Get(); }
	
	/** Set the lower value of the interval. */
	virtual void SetLow( Self* node ) { TreeNode< ValueType >::m_Data =  node->Get(); m_Low = node; }


	/** Get the upper value of the interval. */
	virtual ValueType GetUp() const { return m_Up->Get(); };
	
	/** Set the upper value of the interval. */
	virtual void SetUp( Self* node ) {  m_Up = node; }

	/** Returns the pointer (singleton) of the leaf associated 
	 * with the lower value of the interval.
	 */
	template< class TSelf = Self >
	TSelf* GetLowLeaf() { return dynamic_cast< TSelf* >( m_Low ); }

	/** Returns the pointer (singleton) of the leaf associated 
	 * with the upper value of the interval.
	 */
	template< class TSelf = Self >
	TSelf* GetUpLeaf() { return dynamic_cast< TSelf* >( m_Up ); }

	/** Right rotation. This must be called from a root node, and it
	 * returns the root, so it can be updated if it has been changed.
	 * The node parameter is the location where we do the rotation.
	 */
	template< class TSelf = Self >
	TSelf* RotateRight( TSelf* node );

	/** Right rotation. This must be called from a root node, and it
	 * returns the root, so it can be updated if it has been changed.
	 * The node parameter is the location where we do the rotation.
	 */
	template< class TSelf = Self >
	TSelf* RotateLeft( TSelf* node );

	/** Equilibrate the tree from a node. Should be called from the root.
	 * It returns the new root in case it changed. The recursive parameter
	 * tells if we equilibrate the whole tree, or just the first unbalanced
	 * node.
	 */
	template< class TSelf = Self >
	TSelf* Equilibrate( TSelf *node, bool recursive = false );

	/** Returns the child of the node. The parameter should be Right or Left. */
	template< class TSelf = Self >
	TSelf* GetChild( ChildIdentifier number ) const { return Superclass::template GetChild< TSelf >( number ); }

	/** Returns the parent of the node. If the node is a root, it returns a NULL pointer. */
	template< class TSelf = Self >
	TSelf* GetParent() const { return Superclass::template GetParent< TSelf >(); }

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


	/** Returns true if self object's upper bound is lower than node's lower bound . */
	template< class TSelf = Self >
	bool IsOutsideUpper( TSelf* node ) const;

	/** Returns true if self object's lower bound is upper than node's upper bound . */
	template< class TSelf = Self >
	bool IsOutsideLower( TSelf* node ) const;

	/** Returns right pointer's node to element node.
	 * It should be used to go from a leaf to another.
	 */
	template< class TSelf = Self >
	TSelf* GetRight( const ValueType element );

	/** Returns right pointer's node of itself. It should be used to go from a leaf to another. */
	template< class TSelf = Self >
	TSelf* GetRight();

	/** Returns left pointer's node to element node.
	 * It should be used to go from a leaf to another.
	 */
	template< class TSelf = Self >
	TSelf* GetLeft( const ValueType element );

	/** Returns left pointer's node of itself. It should be used to go from a leaf to another. */
	template< class TSelf = Self >
	TSelf* GetLeft();

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
	TSelf* Split( TSelf *node, TSelf **root );

	/** Copies other's interval to itself. */
	virtual void Set( Self *other );

	/** Sets element to itself. Should only be used for singletons. */
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
	AVLIntervalTreeNode();
	~AVLIntervalTreeNode(){};

	/** Pointers to the upper and lower leaf extremities of the interval. */
	 Self *m_Up, *m_Low;

private:
	AVLIntervalTreeNode( const Self & ); // purposely node implemented
	void operator=( const Self & ); // purposely node implemented

};

template< class T >
std::ostream& operator<<( std::ostream& flux, AVLIntervalTreeNode< T > const& node );

}// namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAVLIntervalTreeNode.hxx"
#endif

#endif
