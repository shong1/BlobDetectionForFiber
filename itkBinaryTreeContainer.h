#ifndef __itkBinaryTreeContainer_h
#define __itkBinaryTreeContainer_h

#include "itkTreeContainerBase.h"
#include "itkTreeNode.h"

namespace itk
{
/** \class BinaryTreeContainer
 * \brief Base class for binary tree manipulations.
 *
 * This class creates a binary tree from some given data.
 * The template determines the type of considered binary tree.
 *
 */
template< class TTreeNode >
class BinaryTreeContainer : public TreeContainerBase< typename TTreeNode::ValueType >
{
public:
	/** Standard typedefs. */
	typedef BinaryTreeContainer									Self;
	typedef SmartPointer< Self >								Pointer;
	typedef SmartPointer< const Self >							ConstPointer;
	typedef TreeContainerBase< typename TTreeNode::ValueType >	Superclass;

	/** Method for creation through the object factory. */
	itkNewMacro( Self );

	/** Run-time type information (and related methods). */
	itkTypeMacro( BinaryTreeContainer, TreeContainerBase );

	/** Convenient typedefs. */
	typedef TTreeNode								TreeNodeType;
	typedef typename TreeNodeType::ValueType		ValueType;
	typedef typename TreeNodeType::Pointer			TreeNodePointer;
	typedef typename TreeNodeType::ChildIdentifier	ChildIdentifier;
	typedef TreeNode< ValueType >					TreeNodeBaseType;

	/** Right and Left indexes. */
	static const ChildIdentifier Right = TreeNodeType::Right, Left = TreeNodeType::Left;

	/** Gets rid of evey nodes of the tree. */
	virtual bool Clear();
	virtual bool Clear( TreeNodeBaseType * node );

	/** Tells if the tree contains an element. */
	virtual bool Contains( const ValueType element );

	/** Counts the number of nodes in the tree. */
	virtual int Count() const;
	virtual int Count( TreeNodeType * node ) const;

	/** Returns the root of the tree. */
	TreeNodeBaseType * GetRoot() const;
	
	/** Tells if the element is a leaf. */
	virtual bool IsLeaf( const ValueType element );

	/** Tells if the root is element. */
	virtual bool IsRoot( const ValueType element );

	/** Sets the root of the tree. */
	virtual bool SetRoot( const ValueType element );
	virtual bool SetRoot( TreeNodeBaseType * node );

	/** Add an element to the tree. */
	virtual TreeNodeType* Add( const ValueType element );

	/** Add an allocated node to the tree. */
	virtual void Add( TreeNodeType *node );

	/** Remove an element from the tree. */
	virtual TreeNodeType* Remove( const ValueType element );

	/** Remove a node with known pointer of the tree. */
	virtual void Remove( TreeNodeType *node );

	/** Splits the tree in two parts. The element will be in the right side.
	 * The returned node is the left side. 
	 */
	virtual TreeNodeType* Split( const ValueType element );
	
	/** Splits the tree in two parts. The pointer node is in the right side. 
	 * The returned node is the left side. 
	 */
	virtual TreeNodeType* Split( TreeNodeType *node );

	/** Splits the tree in two parts, and puts the left output in the allocated tree.
	 * The pointer node is in the right side, tree is the left side.
	 */
	virtual TreeNodeType* Split( TreeNodeType *node, Self *tree );
	
	/** Splits the tree in two parts, and puts the left output in the root pointer.
	 * The pointer node is in the right side.
	 */
	virtual TreeNodeType* Split( TreeNodeType *node, TreeNodeType **root );

	/** Splits the tree in two parts, and puts the left output in the root pointer.
	 * The element is in the right side.
	 */
	virtual TreeNodeType* Split( const ValueType element, Self *tree );

	/** Merges the tree of root node with the self object. The self object
	 * is modified.
	 */
	virtual void Merge( TreeNodeType *node );

	/** Merges the tree with the self object. The tree becomes empty and the self
	 * object is the merged tree.
	 */
	virtual void Merge( Self *tree );


protected:
	BinaryTreeContainer(){ m_Root = NULL; };
	virtual ~BinaryTreeContainer(){};

	TreeNodePointer m_Root;

private:
	BinaryTreeContainer( const Self & );	// purposely not implemented
	void operator=( const Self & );			// purposely not implemented

};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBinaryTreeContainer.hxx"
#endif

#endif
