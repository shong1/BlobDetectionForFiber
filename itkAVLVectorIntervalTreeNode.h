#ifndef __itkAVLVectorIntervalTreeNode_h
#define __itkAVLVectorIntervalTreeNode_h

#include "itkAVLIntervalTreeNode.h"

namespace itk
{
/** \class AVLVectorIntervalTreeNode
 * \brief AVL vector interval tree node implementation, to be used with an herited class
 * from BinaryTreeContainer
 *
 * All the behavior of this class is similar with AVLIntervalTreeNode class.
 * The only difference is the possibility to use a vector, and then use one
 * dimension of it.
 *
 */

template< class TValue >
class AVLVectorIntervalTreeNode : public AVLIntervalTreeNode< TValue >
{
public:
	typedef AVLVectorIntervalTreeNode				Self;
	typedef SmartPointer< Self >					Pointer;
	typedef SmartPointer< const Self >				ConstPointer;
	typedef AVLIntervalTreeNode< TValue >			Superclass;

	itkNewMacro( Self );

	itkTypeMacro( AVLVectorIntervalTreeNode, AVLIntervalTreeNode );

	typedef TValue									ValueType;

	typedef typename Superclass::ChildIdentifier	ChildIdentifier;
	typedef typename Superclass::ChildrenListType	ChildrenListType;
	typedef typename Superclass::Baseclass			Baseclass;

	static const ChildIdentifier Left = Superclass::Left, Right = Superclass::Right;
	static const int MaxDepth = Superclass::MaxDepth;

	unsigned int* GetDimension(){ return m_Dimension; }
	void SetDimension( unsigned int *dimension ){ m_Dimension = dimension; }
	unsigned int* GetOtherDimension(){ return m_OtherDimension; }
	void SetOtherDimension( unsigned int *dimension ){ m_OtherDimension = dimension; }

	/** Get the duplicate node in case there are many nodes with the same element. */
	Self* GetDuplicate(){ return m_Duplicate; }
	/** Set the duplicate node in case there are many nodes with the same element. */
	void SetDuplicate( Self *duplicate ){ m_Duplicate = duplicate; }

	template< class TSelf = Self >
	TSelf* GetLowLeaf() { return Superclass::template GetLowLeaf< TSelf >(); }
	template< class TSelf = Self >
	TSelf* GetUppLeaf() { return Superclass::template GetUppLeaf< TSelf >(); }

	template< class TSelf = Self >
	TSelf* RotateRight( TSelf* node ){ return Superclass::template RotateRight< TSelf >( node ); }
	template< class TSelf = Self >
	TSelf* RotateLeft( TSelf* node ){ return Superclass::template RotateLeft< TSelf >( node ); }
	template< class TSelf = Self >
	TSelf* Equilibrate( TSelf *node, bool recursive = false ){ return Superclass::template Equilibrate< TSelf >( node, recursive ); }

	template< class TSelf = Self >
	TSelf* GetChild( ChildIdentifier number ) const { return Superclass::template GetChild< TSelf >( number ); }
	template< class TSelf = Self >
	TSelf* GetParent() const { return Superclass::template GetParent< TSelf >(); }

	template< class TSelf = Self >
	TSelf* GetNode( const ValueType element ){ return Superclass::template GetNode< TSelf >( element ); }
	template< class TSelf = Self >
	TSelf* GetLowerNode( const ValueType element ){ return Superclass::template GetLowerNode< TSelf >( element ); }
	template< class TSelf = Self >
	TSelf* GetUpperNode( const ValueType element ){ return Superclass::template GetUpperNode< TSelf >( element ); }
	template< class TSelf = Self >
	TSelf* GetStrictlyLowerNode( const ValueType element ){ return Superclass::template GetStrictlyLowerNode< TSelf >( element ); }
	template< class TSelf = Self >
	TSelf* GetStrictlyUpperNode( const ValueType element ){ return Superclass::template GetStrictlyUpperNode< TSelf >( element ); }
	template< class TSelf = Self >
	TSelf* GetLowestNode( int depth = MaxDepth ){ return Superclass::template GetLowestNode< TSelf >( depth ); }
	template< class TSelf = Self >
	TSelf* GetUppestNode( int depth = MaxDepth ){ return Superclass::template GetUppestNode< TSelf >( depth ); }
	template< class TSelf = Self >
	TSelf* GetLowerLeaf( const ValueType element ){ return Superclass::template GetLowerLeaf< TSelf >( element ); }
	template< class TSelf = Self >
	TSelf* GetUpperLeaf( const ValueType element ){ return Superclass::template GetUpperLeaf< TSelf >( element ); }
	template< class TSelf = Self >
	TSelf* GetStrictlyLowerLeaf( const ValueType element ){ return Superclass::template GetStrictlyLowerLeaf< TSelf >( element ); }
	template< class TSelf = Self >
	TSelf* GetStrictlyUpperLeaf( const ValueType element ){ return Superclass::template GetStrictlyUpperLeaf< TSelf >( element ); }

	template< class TSelf = Self >
	bool IsOutsideUpper( TSelf* node ) const;
	template< class TSelf = Self >
	bool IsOutsideLower( TSelf* node ) const;

	/** Duplicate the node and sorts the array of duplicated in function
	 * of the OtherDimension attribute.
	 */
	template< class TSelf = Self >
	TSelf* Duplicate( TSelf *node, ChildIdentifier side );

	template< class TSelf = Self >
	TSelf* GetRight( const ValueType element ){ return Superclass::template GetRight< TSelf >( element ); }
	template< class TSelf = Self >
	TSelf* GetRight(){ return Superclass::template GetRight< TSelf >(); }
	template< class TSelf = Self >
	TSelf* GetLeft( const ValueType element ){ return Superclass::template GetLeft< TSelf >( element ); }
	template< class TSelf = Self >
	TSelf* GetLeft(){ return Superclass::template GetLeft< TSelf >(); }

	virtual void ShareInfo( Self* node ){ m_OtherDimension = node->GetOtherDimension(); m_Dimension = node->GetDimension(); }

	template< class TSelf = Self >
	TSelf* Add( Self *node ){ return Superclass::template Add< TSelf >( node ); }
	template< class TSelf = Self >
	TSelf* Remove( Self *node ){ return Superclass::template Remove< TSelf >( node ); }
	template< class TSelf = Self >
	TSelf* Merge( Self *root ){ return Superclass::template Merge< TSelf >( root ); }
	template< class TSelf = Self >
	TSelf* Split( Self *node, TSelf **root ){ return Superclass::template Split< TSelf >( node, root ); }

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
	AVLVectorIntervalTreeNode(){ m_Duplicate = NULL; }
	~AVLVectorIntervalTreeNode(){};

	unsigned int *m_Dimension, *m_OtherDimension;

	Self *m_Duplicate;


private:
	AVLVectorIntervalTreeNode( const Self & ); // purposely node implemented
	void operator=( const Self & ); // purposely node implemented

};

}// namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAVLVectorIntervalTreeNode.hxx"
#endif

#endif
