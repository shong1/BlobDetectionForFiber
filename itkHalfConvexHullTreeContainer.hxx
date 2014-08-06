#ifndef __itkHalfConvexHullTreeContainer_hxx
#define __itkHalfConvexHullTreeContainer_hxx

#include "itkHalfConvexHullTreeContainer.h"

namespace itk
{

template< class TTreeNode >
typename HalfConvexHullTreeContainer< TTreeNode >
::ConcatenableQueueNodeType* HalfConvexHullTreeContainer< TTreeNode >
::GetQRoot()
{
	return dynamic_cast< TreeNodeType* >( this->GetRoot() )->GetQRoot();
}

template< class TTreeNode >
bool HalfConvexHullTreeContainer< TTreeNode >
::Clear()
{
	return this->Clear( m_Root );
}

template< class TTreeNode >
bool HalfConvexHullTreeContainer< TTreeNode >
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
bool HalfConvexHullTreeContainer< TTreeNode >
::Contains( const ValueType element )
{
	TreeNodeType *node = dynamic_cast< TreeNodeType* >( this->GetRoot() );
	while( !( node->GetQ()->Contains( element ) ) )
	{
		if( node->GetChild( Left ) != NULL && node->GetChild( Left )->GetQRoot()->GetUpp()[m_Dimension] >= element[m_Dimension] )
		{
			node = node->Down( node->GetChild( Left ) );
		}
		else if( node->GetChild( Right ) != NULL && node->GetChild( Right )->GetQRoot()->GetLow()[m_Dimension] <= element[m_Dimension] )
		{
			node = node->Down( node->GetChild( Right ) );
		}
		else
		{
			dynamic_cast< TreeNodeType* >( this->GetRoot() )->Up( node );
			return false;
		}
	}
	dynamic_cast< TreeNodeType* >( this->GetRoot() )->Up( node );
	return true;
}

template< class TTreeNode >
int HalfConvexHullTreeContainer< TTreeNode >
::Count() const
{
	return this->Count( m_Root );
}

template< class TTreeNode >
int HalfConvexHullTreeContainer< TTreeNode >
::Count( TreeNodeType * node ) const
{
	if( node == NULL )
	{
		return 0;
	}
	return node->GetSize();
}

template< class TTreeNode >
typename HalfConvexHullTreeContainer< TTreeNode >
::TreeNodeBaseType * HalfConvexHullTreeContainer< TTreeNode >
::GetRoot() const
{
	return dynamic_cast< TreeNodeBaseType* >( m_Root.GetPointer() );
}

template< class TTreeNode >
bool HalfConvexHullTreeContainer< TTreeNode >
::IsLeaf( const ValueType element )
{
	return true;
}

template< class TTreeNode >
bool HalfConvexHullTreeContainer< TTreeNode >
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
bool HalfConvexHullTreeContainer< TTreeNode >
::SetRoot( const ValueType element )
{
	m_Root = TreeNodeType::New();
	m_Root.GetPointer()->Register();
	m_Root->Set( element );
	return true;
}

template< class TTreeNode >
bool HalfConvexHullTreeContainer< TTreeNode >
::SetRoot( TreeNodeBaseType * node )
{
	m_Root = dynamic_cast< TreeNodeType * >( node );
	if( m_Root.IsNull() )
	{
		return true;
	}
	return true;
}

template< class TTreeNode >
typename HalfConvexHullTreeContainer< TTreeNode >
::ConcatenableQueueNodeType* HalfConvexHullTreeContainer< TTreeNode >
::Add( const VectorType element )
{
	typename ConcatenableQueueNodeType::Pointer node = ConcatenableQueueNodeType::New();
	node.GetPointer()->Register();
	node->Set( element );
	bool newRoot = ( this->GetRoot() == NULL );
	this->Add( node.GetPointer() );
	std::cout << " added " << std::endl;
	if( newRoot )
	{
		dynamic_cast< TreeNodeType * >( this->GetRoot() )->SetTree( this );
	}
	return node;
}

template< class TTreeNode >
void HalfConvexHullTreeContainer< TTreeNode >
::Add( ConcatenableQueueNodeType *node )
{
	node->SetDimension( &m_Dimension );
	node->SetOtherDimension( &m_OtherDimension );
	this->SetRoot( m_Root->Add( node ) );
}

template< class TTreeNode >
typename HalfConvexHullTreeContainer< TTreeNode >
::ConcatenableQueueNodeType* HalfConvexHullTreeContainer< TTreeNode >
::Remove( const VectorType element )
{
	if( this->GetRoot() != NULL )
	{
		this->SetRoot( dynamic_cast< TreeNodeType* >( this->GetRoot() )->Remove( element ) );
	}
	return NULL;
}


}// namespace itk

#endif
