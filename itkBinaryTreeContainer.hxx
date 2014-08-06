#ifndef __itkBinaryTreeContainer_hxx
#define __itkBinaryTreeContainer_hxx

#include "itkBinaryTreeContainer.h"

namespace itk
{

template< class TTreeNode >
bool BinaryTreeContainer< TTreeNode >
::Clear()
{
	return this->Clear( m_Root );
}

template< class TTreeNode >
bool BinaryTreeContainer< TTreeNode >
::Clear( TreeNodeBaseType * node )
{
	if( node == NULL )
	{
		return true;
	}
	if( node->HasChildren() )
	{
		bool success = this->Clear( node->GetChild( Left ) ) && this->Clear( node->GetChild( Right ) );
		success &= node->Remove( node->GetChild( Left ) );
		success &= node->Remove( node->GetChild( Left ) );
		return success;
	}
	return true;
}

template< class TTreeNode >
bool BinaryTreeContainer< TTreeNode >
::Contains( const ValueType element )
{
	return m_Root->template GetNode< TreeNodeType >( element ) != NULL;
}

template< class TTreeNode >
int BinaryTreeContainer< TTreeNode >
::Count() const
{
	return this->Count( m_Root );
}

template< class TTreeNode >
int BinaryTreeContainer< TTreeNode >
::Count( TreeNodeType * node ) const
{
	if( node == NULL )
	{
		return 0;
	}
	return node->GetSize();
}

template< class TTreeNode >
typename BinaryTreeContainer< TTreeNode >
::TreeNodeBaseType * BinaryTreeContainer< TTreeNode >
::GetRoot() const
{
	return dynamic_cast< TreeNodeBaseType* >( m_Root.GetPointer() );
}

template< class TTreeNode >
bool BinaryTreeContainer< TTreeNode >
::IsLeaf( const ValueType element )
{
	const TreeNodeType * node = m_Root->template GetNode< TreeNodeType >( element );
	if( node == NULL )
	{
		return false;
	}
	return !node->HasChildren();
}

template< class TTreeNode >
bool BinaryTreeContainer< TTreeNode >
::IsRoot( const ValueType element )
{
	if( m_Root.IsNull() )
	{
		return false;
	}
	if( m_Root->Get() == element )
	{
		return true;
	}
}

template< class TTreeNode >
bool BinaryTreeContainer< TTreeNode >
::SetRoot( const ValueType element )
{
	m_Root = TreeNodeType::New();
	m_Root.GetPointer()->Register();
	m_Root->Set( element );
	return true;
}

template< class TTreeNode >
bool BinaryTreeContainer< TTreeNode >
::SetRoot( TreeNodeBaseType * node )
{
	std::cout << " SetRoot " << std::endl;
	if( node == NULL )
	{
		std::cout << " ici " << std::endl;
		m_Root = NULL;
	}
	if( m_Root.GetPointer() == node )
	{
		return true;
	}
	m_Root = dynamic_cast< TreeNodeType* >( node );
	if( m_Root.IsNull() )
	{
		return true;
	}
	return true;
}

template< class TTreeNode >
typename BinaryTreeContainer< TTreeNode >
::TreeNodeType * BinaryTreeContainer< TTreeNode >
::Add( const ValueType element )
{
	TreeNodePointer node = TreeNodeType::New();
	node.GetPointer()->Register();
	node->Set( element );
	this->Add( node.GetPointer() );
	if( node != this->GetRoot() )
	{
		node->ShareInfo( m_Root );
	}
	return node.GetPointer();
}

template< class TTreeNode >
void BinaryTreeContainer< TTreeNode >
::Add( TreeNodeType *node )
{
	this->SetRoot( m_Root->Add( node ) );
}

template< class TTreeNode >
typename BinaryTreeContainer< TTreeNode >
::TreeNodeType* BinaryTreeContainer< TTreeNode >
::Remove( const ValueType element )
{
	TreeNodeType *node = m_Root->template GetNode< TreeNodeType >( element );
	std::cout << node << std::endl;
	this->Remove( node );
	return node;
}

template< class TTreeNode >
void BinaryTreeContainer< TTreeNode >
::Remove( TreeNodeType *node )
{
	this->SetRoot( m_Root->Remove( node ) );
}

template< class TTreeNode >
typename BinaryTreeContainer< TTreeNode >
::TreeNodeType* BinaryTreeContainer< TTreeNode >
::Split( const ValueType element )
{
	return this->Split( m_Root->GetNode( element ) );
}

template< class TTreeNode >
typename BinaryTreeContainer< TTreeNode >
::TreeNodeType* BinaryTreeContainer< TTreeNode >
::Split( TreeNodeType *node )
{
	TreeNodeType *root = NULL;
	return m_Root->Split( node, &root );
}

template< class TTreeNode >
typename BinaryTreeContainer< TTreeNode >
::TreeNodeType* BinaryTreeContainer< TTreeNode >
::Split( const ValueType element, Self* tree )
{
	return this->Split( m_Root->GetNode( element ), tree );
}

template< class TTreeNode >
typename BinaryTreeContainer< TTreeNode >
::TreeNodeType* BinaryTreeContainer< TTreeNode >
::Split( TreeNodeType *node, Self* tree )
{
	std::cout << " this split " << std::endl;
	if( tree == NULL )
	{
		return m_Root.GetPointer();
	}
	TreeNodeType *root = dynamic_cast< TreeNodeType* >( tree->GetRoot() );
	this->SetRoot( m_Root->Split( node, &root ) );
	tree->SetRoot( root );
	return root;
}

template< class TTreeNode >
typename BinaryTreeContainer< TTreeNode >
::TreeNodeType* BinaryTreeContainer< TTreeNode >
::Split( TreeNodeType *node, TreeNodeType **root )
{
	this->SetRoot( m_Root->Split( node, root ) );
	return m_Root.GetPointer();
}

template< class TTreeNode >
void BinaryTreeContainer< TTreeNode >
::Merge( TreeNodeType *node )
{
	this->SetRoot( m_Root->Merge( node ) );
}

template< class TTreeNode >
void BinaryTreeContainer< TTreeNode >
::Merge( Self* tree )
{
	if( tree == NULL )
	{
		return;
	}
	this->Merge( dynamic_cast< TreeNodeType* >( tree->GetRoot() ) );
	tree->SetRoot( nullptr );
}

}// namespace itk

#endif
