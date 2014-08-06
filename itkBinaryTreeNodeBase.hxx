#ifndef __itkBinaryTreeNodeBase_hxx
#define __itkBinaryTreeNodeBase_hxx

#include "itkBinaryTreeNodeBase.h"

namespace itk
{

template< class TValue >
BinaryTreeNodeBase< TValue >::BinaryTreeNodeBase() : m_Depth( 1 ), m_Size( 1 )
{}

template< class TValue >
void BinaryTreeNodeBase< TValue >
::Set( const Self & other )
{
	this->Set( other.Get() );
}

template< class TValue >
void BinaryTreeNodeBase< TValue >
::Set( const ValueType & other )
{
	Superclass::m_Data = other;
}

template< class TValue >
void BinaryTreeNodeBase< TValue >
::SetParent( Self * node )
{
	Superclass::m_Parent = node;
}

template< class TValue >
void BinaryTreeNodeBase< TValue >
::AddChild( ChildIdentifier position, Self *node )
{
	if( this == NULL )
	{
		return;
	}
	if( position > 1 )
	{
		std::cout << "Error: BinaryTreeNodeBase< TValue >::AddChild( position, node ) : position > 1" << std::endl;
	}
	if( !Superclass::m_Children.size() )
	{
		if( node != NULL )
		{
			Superclass::m_Children.push_back( NULL ); 
			Superclass::m_Children.push_back( NULL ); 
		}
		else
		{
			return;
		}
	}
	Superclass::m_Children[position] = node;
	if( node != NULL )
	{
		node->SetParent( this );
	}
	if( node == NULL )
	{
		if( Superclass::m_Children[Left].IsNull() && Superclass::m_Children[Right].IsNull() )
		{
			Superclass::m_Children.clear();
		}
	}
}

template< class TValue >
template< class TSelf >
TSelf* BinaryTreeNodeBase< TValue >
::GetNode( const ValueType element ) const
{
	if( this == NULL )
	{
		return NULL;
	}
	if( *this == element )
	{
		return dynamic_cast< TSelf* >( this );
	}
	TSelf* tmp = this->GetChild( Left )->template GetNode< TSelf >( element );
	if( tmp == NULL )
	{
		return this->GetChild( Right )->template GetNode< TSelf >( element );
	}
	return tmp;
}

template< class TValue >
template< class TSelf >
TSelf* BinaryTreeNodeBase< TValue >
::RotateRight( TSelf* node )
{
	if( node == NULL )
	{
		return NULL;
	}
	TSelf *parent = node->template GetParent< TSelf >();
	TSelf *newRoot = dynamic_cast< TSelf* >( this );
	TSelf *pivot = node->template GetChild< TSelf >( Left );
	node->AddChild( Left, pivot->GetChild( Right ) );
	pivot->AddChild( Right, node );
	if( parent == NULL )
	{
		pivot->SetParent( NULL );
		newRoot = pivot;
	}
	else
	{
		parent->AddChild( parent->ChildPosition( node ), pivot );
	}
	node->UpdateDepthToRoot( newRoot );
	return newRoot;
}

template< class TValue >
template< class TSelf >
TSelf* BinaryTreeNodeBase< TValue >
::RotateLeft( TSelf* node )
{
	if( node == NULL )
	{
		return NULL;
	}
	TSelf *parent = node->template GetParent< TSelf >();
	TSelf *newRoot = dynamic_cast< TSelf*>( this );
	TSelf *pivot = node->template GetChild< TSelf >( Right );
	node->AddChild( Right, pivot->GetChild( Left ) );
	pivot->AddChild( Left, node );
	if( parent == NULL )
	{
		pivot->SetParent( NULL );
		newRoot = pivot;
	}
	else
	{
		parent->AddChild( parent->ChildPosition( node ), pivot );
	}
	node->UpdateDepthToRoot( newRoot );
	return newRoot;
}

template< class TValue >
template< class TSelf >
TSelf* BinaryTreeNodeBase< TValue >
::Rotate( TSelf * root, ChildIdentifier side )
{
	switch( side )
	{
		case Right:
			return this->template RotateRight< TSelf >( root );
			break;
		case Left:
			return this->template RotateLeft< TSelf >( root );
			break;
		default:
			return root;
			break;
	}
}

template< class TValue >
void BinaryTreeNodeBase< TValue >
::WipeChildren()
{
	Superclass::m_Children.clear();
}

template< class TValue >
void BinaryTreeNodeBase< TValue >
::UpdateDepthToRoot( Self *root )
{
	if( root == NULL )
    {
        return;
	}
	if( !this->HasChildren() )
	{
		m_Depth = 1;
		m_Size = 1;
		if( this->GetParent() != NULL )
		{
			this->GetParent()->UpdateDepthToRoot( root  );
		}
		return;
	}
	Self *child;
	int dleft = dleft = this->GetChild( Left )->GetDepth();
	int dright = this->GetChild( Right )->GetDepth();
	if( dleft > dright )
	{
		child = this->GetChild( Left );
	}
	else
	{
		child = this->GetChild( Right );
	}
	Self *node = this, *other;
	node->SetDepth( child->GetDepth() + 1 );
	node->SetSize( node->GetChild( Right )->GetSize() + node->GetChild( Left )->GetSize() + 1 );
	child = node;
	node = node->GetParent();
	bool updateDepth = true;
	while( node != NULL  )
	{
		other = node->GetChild( !node->ChildPosition( child ) );
		if( updateDepth )
		{
			if( other->GetDepth() > child->GetDepth() )
			{
				if( node->GetDepth() == other->GetDepth() + 1 )
				{
					updateDepth = false;
				}
				node->SetDepth( other->GetDepth() + 1 );
			}
			else
			{
				node->SetDepth( other->GetDepth() + 1 );
				node->SetDepth( child->GetDepth() + 1 );
			}
		}
		node->SetSize( node->GetChild( Right )->GetSize() + node->GetChild( Left )->GetSize() + 1 );
		child = node;
		node = node->GetParent();
	}
	if( node != NULL )
	{

	}
	if( node == this )
	{
		std::cout << "error: Binary Tree is relooping on itself" << std::endl;
	}
}

template< class TValue >
int BinaryTreeNodeBase< TValue >
::GetSize() const
{
	if( this == NULL )
	{
		return 0;
	}
	return m_Size;
}

template< class TValue >
int BinaryTreeNodeBase< TValue >
::GetDepth() const
{
	if( this == NULL )
	{
		return 0;
	}
	return m_Depth;
}

template< class TValue >
void BinaryTreeNodeBase< TValue >
::IncrementDepth()
{
	++m_Depth;
}

template< class TValue >
void BinaryTreeNodeBase< TValue >
::DecrementDepth()
{
	--m_Depth;
}

template< class TValue >
template< class TSelf >
bool BinaryTreeNodeBase< TValue >
::operator<( const TSelf &other )
{
	return this->Get() < other.Get();
}

template< class TValue >
template< class TSelf >
bool BinaryTreeNodeBase< TValue >
::operator>( const TSelf &other )
{
	return this->Get() > other.Get();
}

template< class TValue >
template< class TSelf >
bool BinaryTreeNodeBase< TValue >
::operator<=( const TSelf &other )
{
	return this->Get() <= other.Get();
}

template< class TValue >
template< class TSelf >
bool BinaryTreeNodeBase< TValue >
::operator>=( const TSelf &other )
{
	return this->Get() >= other.Get();
}

template< class TValue >
template< class TSelf >
bool BinaryTreeNodeBase< TValue >
::operator==( const TSelf &other )
{
	return this->Get() == other.Get();
}

template< class TValue >
template< class TSelf >
bool BinaryTreeNodeBase< TValue >
::operator!=( const TSelf &other )
{
	return this->Get() != other.Get();
}

template< class TValue >
template< class TSelf >
bool BinaryTreeNodeBase< TValue >
::operator<( const ValueType & element )
{
	return this->Get() < element;
}

template< class TValue >
template< class TSelf >
bool BinaryTreeNodeBase< TValue >
::operator>( const ValueType & element )
{
	return this->Get() > element;
}

template< class TValue >
template< class TSelf >
bool BinaryTreeNodeBase< TValue >
::operator<=( const ValueType & element )
{
	return this->Get() <= element;
}

template< class TValue >
template< class TSelf >
bool BinaryTreeNodeBase< TValue >
::operator>=( const ValueType & element )
{
	return this->Get() >= element;
}

template< class TValue >
template< class TSelf >
bool BinaryTreeNodeBase< TValue >
::operator==( const ValueType & element )
{
	return this->Get() == element;
}

template< class TValue >
template< class TSelf >
bool BinaryTreeNodeBase< TValue >
::operator!=( const ValueType & element )
{
	return this->Get() != element;
}

}// namespace itk

#endif
