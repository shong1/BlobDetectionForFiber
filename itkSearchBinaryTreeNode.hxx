#ifndef __itkSearchBinaryTreeNode_hxx
#define __itkSearchBinaryTreeNode_hxx

#include "itkSearchBinaryTreeNode.h"

namespace itk
{

template< class TValue >
SearchBinaryTreeNode< TValue >::SearchBinaryTreeNode()
{
}

template< class TValue >
template< class TSelf >
TSelf* SearchBinaryTreeNode< TValue >
::Add( TSelf *tmp )
{
	std::cout << " here " << std::endl;
	if( tmp == NULL )
	{
		return dynamic_cast< TSelf* >( this );
	}
	if( this == NULL )
	{
		return tmp;
	}
	TSelf *node = dynamic_cast< TSelf* >( this );
	while( 1 )
	{
		if( *node < *tmp )
		{
			if( !node->HasChildren() )
			{
				node->AddChild( Right, tmp );
				node->UpdateDepthToRoot( this );
				return dynamic_cast< TSelf* >( this );
			}
			if( node->GetChild( Right ) == NULL )
			{
				node->AddChild( Right, tmp );
				node->UpdateDepthToRoot( this );
				return dynamic_cast< TSelf* >( this );
			}
			node = node->template GetChild< TSelf >( Right );
		}
		else if( *node > *tmp )
		{
			if( !node->HasChildren() )
			{
				node->AddChild( Left, tmp );
				node->UpdateDepthToRoot( this );
				return dynamic_cast< TSelf* >( this );
			}
			if( node->GetChild( Left ) == NULL )
			{
				node->AddChild( Left, tmp );
				node->UpdateDepthToRoot( this );
				return dynamic_cast< TSelf* >( this );
			}
			node = node->template GetChild< TSelf >( Left );
		}
	}
}

template< class TValue >
template< class TSelf >
TSelf *SearchBinaryTreeNode< TValue >
::Remove( TSelf *node )
{
	if( node == NULL || this == NULL )
	{
		return dynamic_cast< TSelf* >( this );
	}
	Self *parent = node->template GetParent< TSelf >();
	if( !node->HasChildren() )
	{
		if( this == node )
		{
			return NULL;
		}
		parent->AddChild( parent->ChildPosition( node ), NULL );
		parent->UpdateDepthToRoot( this );
		return dynamic_cast< TSelf* >( this );
	}
	if( node->GetChild( Left ) == NULL && node->GetChild( Right ) != NULL )
	{
		if( this == node )
		{
			return node->GetChild( Right );
		}
		parent->AddChild( parent->ChildPosition( node ), node->GetChild( Right ) );
		parent->UpdateDepthToRoot( this );
		return dynamic_cast< TSelf* >( this );
	}
	if( node->GetChild( Left ) != NULL && node->GetChild( Right ) == NULL )
	{
		if( this == node )
		{
			return node->template GetChild< TSelf >( Left );
		}
		parent->AddChild( parent->ChildPosition( node ), node->GetChild( Left ) );
		parent->UpdateDepthToRoot( this );
		return dynamic_cast< TSelf* >( this );
	}
	if( this == node )
	{
		TSelf *tmp = node->GetChild( Left )->GetUppestNode();
		node->Set( *tmp );
		return dynamic_cast< TSelf* >( this )->template Remove< TSelf >( tmp );
	}
	TSelf *R = node, *L = node;
	unsigned int Rn = 0, Ln = 0;
	while( R->HasChildren() && R->GetChild( Left ) != NULL )
	{
		R = R->template GetChild< TSelf >( Left );
		++Rn;
	}
	while( L->HasChildren() && L->GetChild( Right ) != NULL )
	{
		L = L->template GetChild< TSelf >( Right );
		++Ln;
	}
	if( Rn < Ln )
	{
		node->Set( *R );
		parent = R->GetParent();
		parent->AddChild( parent->ChildPosition( R ), NULL );
		parent->UpdateDepthToRoot( this );
		return dynamic_cast< TSelf* >( this );
	}
	node->Set( *L );
	parent = L->template GetParent< TSelf >();
	parent->AddChild( parent->ChildPosition( L ), NULL );
	parent->UpdateDepthToRoot( this );
	return dynamic_cast< TSelf* >( this );
}

template< class TValue >
template< class TSelf >
TSelf * SearchBinaryTreeNode< TValue >
::Split( TSelf *node, TSelf **root )
{
	ValueType element = node->Get();
	TSelf *tmp, *newRoot = NULL;
	node = dynamic_cast< TSelf* >( this );
	*root = NULL;
	while( node != NULL )
	{
		if( *node < element )
		{
			tmp = node->template GetChild< TSelf >( Right );
			*root = (*root)->template Merge< TSelf >( node->template GetChild< TSelf >( Left ) );
			node->WipeChildren();
			node->SetDepth( 1 );
			node->SetSize( 1 );
			std::cout << " Add " << std::endl;
			std::cout << " root " << (*root)->GetLow() << ", " << (*root)->GetUpp() << std::endl;
			*root = (*root)->template Add< TSelf >( node );
			std::cout << " end " << std::endl;
		}
		else
		{
			tmp = node->template GetChild< TSelf >( Left );
			if( newRoot == NULL )
			{
				std::cout << node->GetChild( Right ) << std::endl;
				newRoot = node->template GetChild< TSelf >( Right );
				newRoot->SetParent( NULL );
			}
			else
			{
				newRoot = newRoot->template Merge< TSelf >( node->template GetChild< TSelf >( Right ) );
			}
			node->WipeChildren();
			node->SetDepth( 1 );
			node->SetSize( 1 );
			newRoot = newRoot->template Add< TSelf >( node );
		}
		node = tmp;
	}
	return newRoot;
}

template< class TValue >
template< class TSelf >
TSelf* SearchBinaryTreeNode< TValue >
::Merge( TSelf *root )
{
	if( root == NULL )
	{
		return dynamic_cast< TSelf* >( this );
	}
	if( this == NULL )
	{
		return root;
	}
	TSelf *right = root->template GetChild< TSelf >( Right ), *left = root->template GetChild< TSelf >( Left );
	root->WipeChildren();
	root->SetDepth( 1 );
	root->SetSize( 1 );
	this->template Add< TSelf >( root );
	this->template Merge< TSelf >( right );
	this->template Merge< TSelf >( left );
	return dynamic_cast< TSelf* >( this );
}

template< class TValue >
template< class TSelf >
TSelf* SearchBinaryTreeNode< TValue >
::GetNode( const ValueType element ) 
{
	std::cout << " GetNode " << std::endl;
	if( this == NULL )
	{
		return NULL;
	}
	std::cout << " element " << element << std::endl;
	TSelf* node = dynamic_cast< TSelf* >( this );
	std::cout << " begin " << std::endl;
	while( node->HasChildren() && *node != element )
	{
		std::cout << " node " << *node << std::endl;
		if( *node < element && ( node->GetChild( Left ) == NULL || node->GetChild( Left ) != NULL && *(node->GetChild( Left )) != element ) )
		{
			node = node->GetChild( Right );
		}
		else
		{
			node = node->GetChild( Left );
		}
	}
	std::cout << " end " << std::endl;
	std::cout << " node " << *node << std::endl;
	std::cout << ( *node == element ) << std::endl;
	if( *node == element )
	{
		return node;
	}
	return NULL;
}

template< class TValue >
template< class TSelf >
TSelf* SearchBinaryTreeNode< TValue >
::GetLowerNode( const ValueType element )
{
	if( this == NULL )
	{
		return NULL;
	}
	TSelf *node = dynamic_cast< TSelf* >( this );
	TSelf *lower = node;
	while( node != NULL )
	{
		if( *node > *lower && *node <= element )
		{
			lower = node;
			if( *node == element )
			{
				return lower;
			}
		}
		if( *node > element )
		{
			node = node->GetChild( Left );
		}
		else
		{
			node = node->GetChild( Right );
		}
	}
	return lower;
}

template< class TValue >
template< class TSelf >
TSelf* SearchBinaryTreeNode< TValue >
::GetUpperNode( const ValueType element )
{
	if( this == NULL )
	{
		return NULL;
	}
	TSelf *node = dynamic_cast< TSelf* >( this );
	TSelf *upper = node;
	while( node != NULL )
	{
		std::cout << " upper node " << node->GetLow() << node->GetUpp() << std::endl;
		if( *node < *upper && *node >= element )
		{
			upper = node;
			if( *node == element )
			{
				return upper;
			}
		}
		if( *node <= element )
		{
			node = node->GetChild( Right );
		}
		else
		{
			node = node->GetChild( Left );
		}
	}
	return upper;
}

template< class TValue >
template< class TSelf >
TSelf* SearchBinaryTreeNode< TValue >
::GetStrictlyLowerNode( const ValueType element )
{
	if( this == NULL )
	{
		return NULL;
	}
	TSelf *node = dynamic_cast< TSelf* >( this );
	TSelf *strictlyLower = node;
	while( node != NULL )
	{
		if( *node > *strictlyLower && *node < element )
		{
			strictlyLower = node;
		}
		if( *node >= element )
		{
			node = node->GetChild( Right );
		}
		else
		{
			node = node->GetChild( Left );
		}
	}
	return strictlyLower;
}

template< class TValue >
template< class TSelf >
TSelf* SearchBinaryTreeNode< TValue >
::GetStrictlyUpperNode( const ValueType element )
{
	if( this == NULL )
	{
		return NULL;
	}
	TSelf *node = dynamic_cast< TSelf* >( this );
	TSelf *strictlyUpper = node;
	while( node != NULL )
	{
		if( *node < *strictlyUpper && *node > element )
		{
			strictlyUpper = node;
		}
		if( *node < element )
		{
			node = node->GetChild( Left );
		}
		else
		{
			node = node->GetChild( Right );
		}
	}
	return strictlyUpper;
}

template< class TValue >
template< class TSelf >
TSelf* SearchBinaryTreeNode< TValue >
::GetLowestNode( int depth )
{
	if( this == NULL )
	{
		return NULL;
	}
	TSelf * node = dynamic_cast< TSelf* >( this );
	if( depth == MaxDepth )
	{
		depth = this->GetDepth();
	}
	while( node->GetChild( Left ) != NULL && depth )
	{
		--depth;
		node = node->GetChild( Left );
	}
	return node;
}

template< class TValue >
template< class TSelf >
TSelf* SearchBinaryTreeNode< TValue >
::GetUppestNode( int depth ) 
{
	if( this == NULL )
	{
		return NULL;
	}
	TSelf * node = dynamic_cast< TSelf* >( this );
	if( depth == MaxDepth )
	{
		depth = this->GetDepth();
	}
	while( node->GetChild( Right ) != NULL && depth )
	{
		--depth;
		node = node->GetChild( Right );
	}
	return node;
}

template< class TValue >
template< class TSelf >
TSelf* SearchBinaryTreeNode< TValue >
::GetLowerLeaf( const ValueType element )
{
	if( this == NULL )
	{
		return NULL;
	}
	TSelf *node = dynamic_cast< TSelf* >( this );
	while( node->HasChildren() )
	{
		if( *node > element )
		{
			node = node->GetChild( Left );
		}
		else
		{
			node = node->GetChild( Right );
		}
	}
	return node;
}

template< class TValue >
template< class TSelf >
TSelf* SearchBinaryTreeNode< TValue >
::GetUpperLeaf( const ValueType element )
{
	if( this == NULL )
	{
		return NULL;
	}
	TSelf *node = dynamic_cast< TSelf* >( this );
	std::cout << " element " << element << std::endl;
	while( node->HasChildren() )
	{
		std::cout << " upper node " << node->GetLow() << node->GetUpp() << std::endl;
		std::cout << " left " << node->GetChild( Left ) << " right " << node->GetChild( Right ) << std::endl;
		std::cout << " depth " << node->GetDepth() << std::endl;
		if( node->GetChild( Right ) != NULL ) // ici
		std::cout << " upper node right " << node->GetChild( Right )->GetLow() << node->GetChild( Right )->GetUpp() << std::endl;
		if( node->GetChild( Left ) != NULL ) // ici
		std::cout << " upper node left " << node->GetChild( Left )->GetLow() << node->GetChild( Left )->GetUpp() << std::endl;
		if( *node < element )
		{
			node = node->GetChild( Right );
		}
		else
		{
			node = node->GetChild( Left );
		}
	}
	return node;
}

template< class TValue >
template< class TSelf >
TSelf* SearchBinaryTreeNode< TValue >
::GetStrictlyLowerLeaf( const ValueType element )
{
	if( this == NULL )
	{
		return NULL;
	}
	TSelf *node = dynamic_cast< TSelf* >( this );
	while( node->HasChildren() )
	{
		if( *node >= element )
		{
			node = node->GetChild( Right );
		}
		else
		{
			node = node->GetChild( Left );
		}
	}
	return node;
}

template< class TValue >
template< class TSelf >
TSelf* SearchBinaryTreeNode< TValue >
::GetStrictlyUpperLeaf( const ValueType element )
{
	if( this == NULL )
	{
		return NULL;
	}
	TSelf *node = dynamic_cast< TSelf* >( this );
	while( node->HasChildren() )
	{
		if( *node < element )
		{
			node = node->GetChild( Left );
		}
		else
		{
			node = node->GetChild( Right );
		}
	}
	return node;
}



}// namespace itk

#endif
