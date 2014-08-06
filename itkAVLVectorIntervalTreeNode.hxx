#ifndef __itkAVLVectorIntervalTreeNode_hxx
#define __itkAVLVectorIntervalTreeNode_hxx

#include "itkAVLVectorIntervalTreeNode.h"

namespace itk
{

template< class TValue >
template< class TSelf >
TSelf * AVLVectorIntervalTreeNode< TValue >
::Duplicate( TSelf* node, ChildIdentifier side )
{
	if( this == NULL || *this != *node )
	{
		return NULL;
	}
	TSelf* current = dynamic_cast< TSelf* >( this );
	while( current->Get()[*m_OtherDimension] > node->Get()[*m_OtherDimension] && side == Right
			|| current->Get()[*m_OtherDimension] < node->Get()[*m_OtherDimension] && side == Left )
	{
		if( current->GetDuplicate() == NULL )
		{
			current->SetDuplicate( node );
			return dynamic_cast< TSelf* >( this );
		}
		current = current->GetDuplicate();
	}
	if( current == this )
	{
		node->SetDuplicate( current );
		node->SetParent( current->GetParent() );
		if( node->GetParent() != NULL )
		{
			current->Register();
			node->GetParent()->AddChild( node->GetParent()->ChildPosition( current ), node );
		}
		return node;
	}
	node->SetDuplicate( current->GetDuplicate() );
	current->SetDuplicate( node );
	return dynamic_cast< TSelf* >( this );
}

template< class TValue >
template< class TSelf >
bool AVLVectorIntervalTreeNode< TValue >
::operator<( const TSelf &other )
{
	if( this->GetChild( Left ) != NULL )
	{
		return this->GetChild( Left )->GetUp()[*m_Dimension] <= other.GetUp()[*m_Dimension];
	}
	if( this->GetChild( Right ) != NULL )
	{
		return this->GetChild( Right )->GetLow()[*m_Dimension] > other.GetLow()[*m_Dimension];
	}
	return this->Get()[*m_Dimension] < other.Get()[*m_Dimension];
}

template< class TValue >
template< class TSelf >
bool AVLVectorIntervalTreeNode< TValue >
::operator<=( const TSelf &other )
{
	if( this->GetChild( Left ) != NULL )
	{
		return this->GetChild( Left )->GetUp()[*m_Dimension] <= other.GetUp()[*m_Dimension]
			&& ( this->GetChild( Left )->GetUp()[*m_Dimension] == other.GetUp()[*m_Dimension]
			&& this->GetChild( Left )->GetUp()[*m_OtherDimension] <= other.GetUp()[*m_OtherDimension]
			|| this->GetChild( Left )->GetUp()[*m_Dimension] != other.GetUp()[*m_Dimension] );
	}
	if( this->GetChild( Right ) != NULL )
	{
		return this->GetChild( Right )->GetLow()[*m_Dimension] >= other.GetLow()[*m_Dimension]
			&& ( this->GetChild( Right )->GetLow()[*m_Dimension] == other.GetLow()[*m_Dimension]
			&& this->GetChild( Right )->GetLow()[*m_OtherDimension] >= other.GetLow()[*m_OtherDimension]
			|| this->GetChild( Right )->GetLow()[*m_Dimension] != other.GetLow()[*m_Dimension] );
	}
	return this->Get()[*m_Dimension] <= other.Get()[*m_Dimension]
		&& ( this->Get()[*m_Dimension] == other.Get()[*m_Dimension]
		&& this->Get()[*m_OtherDimension] <= other.Get()[*m_OtherDimension]
		|| this->Get()[*m_Dimension] != other.Get()[*m_Dimension] );
}

template< class TValue >
template< class TSelf >
bool AVLVectorIntervalTreeNode< TValue >
::operator>( const TSelf &other )
{
	if( this->GetChild( Right ) != NULL )
	{
		return this->GetChild( Right )->GetLow()[*m_Dimension] <= other.GetLow()[*m_Dimension]
			&& ( this->GetChild( Right )->GetLow()[*m_Dimension] == other.GetLow()[*m_Dimension]
			&& this->GetChild( Right )->GetLow()[*m_OtherDimension] <= other.GetLow()[*m_OtherDimension]
			|| this->GetChild( Right )->GetLow()[*m_Dimension] != other.GetLow()[*m_Dimension] );
	}
	if( this->GetChild( Left ) != NULL )
	{
		return this->GetChild( Left )->GetUp()[*m_Dimension] > other.GetUp()[*m_Dimension];
	}
	return this->Get()[*m_Dimension] > other.Get()[*m_Dimension];
}

template< class TValue >
template< class TSelf >
bool AVLVectorIntervalTreeNode< TValue >
::operator>=( const TSelf &other )
{
	if( this->GetChild( Right ) != NULL )
	{
		return this->GetChild( Right )->GetLow()[*m_Dimension] <= other.GetLow()[*m_Dimension]
		&& ( this->GetChild( Right )->GetLow()[*m_Dimension] == other.GetLow()[*m_Dimension]
		&& this->GetChild( Right )->GetLow()[*m_OtherDimension] <= other.GetLow()[*m_OtherDimension]
		|| this->GetChild( Right )->GetLow()[*m_Dimension] != other.GetLow()[*m_Dimension] );
	}
	if( this->GetChild( Left ) != NULL )
	{
		return this->GetChild( Left )->GetUp()[*m_Dimension] >= other.GetUp()[*m_Dimension]
		 	&& ( this->GetChild( Left )->GetUp()[*m_Dimension] == other.GetUp()[*m_Dimension]
		 	&& this->GetChild( Left )->GetUp()[*m_OtherDimension] >= other.GetUp()[*m_OtherDimension]
		 	|| this->GetChild( Left )->GetUp()[*m_Dimension] != other.GetUp()[*m_Dimension] );
	}
	return this->Get()[*m_Dimension] >= other.Get()[*m_Dimension]
		&& ( this->Get()[*m_Dimension] == other.Get()[*m_Dimension]
		&& this->Get()[*m_OtherDimension] >= other.Get()[*m_OtherDimension]
		|| this->Get()[*m_Dimension] != other.Get()[*m_Dimension] );
}

template< class TValue >
template< class TSelf >
bool AVLVectorIntervalTreeNode< TValue >
::operator==( const TSelf &other )
{
	return this->GetUp()[*m_Dimension] == other.GetUp()[*m_Dimension] && this->GetLow()[*m_Dimension] == other.GetLow()[*m_Dimension];
}

template< class TValue >
template< class TSelf >
bool AVLVectorIntervalTreeNode< TValue >
::operator!=( const TSelf &other )
{
	return this->GetUp()[*m_Dimension] != other.GetUp()[*m_Dimension] || this->GetLow()[*m_Dimension] != other.GetLow()[*m_Dimension];
}

template< class TValue >
template< class TSelf >
bool AVLVectorIntervalTreeNode< TValue >
::operator<( const ValueType & element )
{
	if( this->GetChild( Left ) != NULL )
	{
		return this->GetChild( Left )->GetUp()[*m_Dimension] <= element[*m_Dimension];
	}
	if( this->GetChild( Right ) != NULL )
	{
		return this->GetChild( Right )->GetLow()[*m_Dimension] > element[*m_Dimension];
	}
	return this->Get()[*m_Dimension] < element[*m_Dimension];
}

template< class TValue >
template< class TSelf >
bool AVLVectorIntervalTreeNode< TValue >
::operator<=( const ValueType & element )
{
	if( this->GetChild( Left ) != NULL )
	{
		return this->GetChild( Left )->GetUp()[*m_Dimension] <= element[*m_Dimension]
			&& ( this->GetChild( Left )->GetUp()[*m_Dimension] == element[*m_Dimension] 
			&& this->GetChild( Left )->GetUp()[*m_OtherDimension] <= element[*m_OtherDimension]
			|| this->GetChild( Left )->GetUp()[*m_Dimension] != element[*m_Dimension] );
	}
	if( this->GetChild( Right ) != NULL )
	{
		return this->GetChild( Right )->GetLow()[*m_Dimension] >= element[*m_Dimension]
			&&( this->GetChild( Right )->GetLow()[*m_Dimension] == element[*m_Dimension]
			&& this->GetChild( Right )->GetLow()[*m_OtherDimension] >= element[*m_OtherDimension]
			|| this->GetChild( Right )->GetLow()[*m_Dimension] != element[*m_Dimension] );
	}
	return this->Get()[*m_Dimension] <= element[*m_Dimension]
		&& ( this->Get()[*m_Dimension] == element[*m_Dimension]
		&& this->Get()[*m_OtherDimension] <= element[*m_OtherDimension]
		|| this->Get()[*m_Dimension] != element[*m_Dimension] );
}

template< class TValue >
template< class TSelf >
bool AVLVectorIntervalTreeNode< TValue >
::operator>( const ValueType & element )
{
	if( this->GetChild( Right ) != NULL )
	{
		return this->GetChild( Right )->GetLow()[*m_Dimension] <= element[*m_Dimension]
			&& ( this->GetChild( Right )->GetLow()[*m_Dimension] == element[*m_Dimension]
			&& this->GetChild( Right )->GetLow()[*m_OtherDimension] <= element[*m_OtherDimension]
			|| this->GetChild( Right )->GetLow()[*m_Dimension] != element[*m_Dimension] );
	}
	if( this->GetChild( Left ) != NULL )
	{
		return this->GetChild( Left )->GetUp()[*m_Dimension] > element[*m_Dimension];
	}
	return this->Get()[*m_Dimension] > element[*m_Dimension];
}

template< class TValue >
template< class TSelf >
bool AVLVectorIntervalTreeNode< TValue >
::operator>=( const ValueType & element )
{
	if( this->GetChild( Right ) != NULL )
	{
		return this->GetChild( Right )->GetLow()[*m_Dimension] <= element[*m_Dimension]
			&& ( this->GetChild( Right )->GetLow()[*m_Dimension] == element[*m_Dimension]
			&& this->GetChild( Right )->GetLow()[*m_OtherDimension] <= element[*m_OtherDimension]
			|| this->GetChild( Right )->GetLow()[*m_Dimension] != element[*m_Dimension] );
	}
	if( this->GetChild( Left ) != NULL )
	{
		return this->GetChild( Left )->GetUp()[*m_Dimension] >= element[*m_Dimension]
			&& ( this->GetChild( Left )->GetUp()[*m_Dimension] == element[*m_Dimension]
			&& this->GetChild( Left )->GetUp()[*m_OtherDimension] >= element[*m_OtherDimension]
			|| this->GetChild( Left )->GetUp()[*m_Dimension] != element[*m_Dimension] );
	}
	return this->Get()[*m_Dimension] >= element[*m_Dimension]
		&& ( this->Get()[*m_Dimension] == element[*m_Dimension]
		&& this->Get()[*m_OtherDimension] >= element[*m_OtherDimension]
		|| this->Get()[*m_Dimension] != element[*m_Dimension] );
}

template< class TValue >
template< class TSelf >
bool AVLVectorIntervalTreeNode< TValue >
::operator==( const ValueType & element )
{
	return this->GetUp()[*m_Dimension] == element[*m_Dimension] && this->GetLow()[*m_Dimension] == element[*m_Dimension];
}

template< class TValue >
template< class TSelf >
bool AVLVectorIntervalTreeNode< TValue >
::operator!=( const ValueType & element )
{
	return this->GetUp()[*m_Dimension] != element[*m_Dimension] || this->GetLow()[*m_Dimension] != element[*m_Dimension];
}


template< class TValue >
template< class TSelf >
bool AVLVectorIntervalTreeNode< TValue >
::IsOutsideUpper( TSelf* node ) const
{
	return this->GetUp()[*m_Dimension] < node->GetLow()[*m_Dimension];
}

template< class TValue >
template< class TSelf >
bool AVLVectorIntervalTreeNode< TValue >
::IsOutsideLower( TSelf* node ) const
{
	return this->GetLow()[*m_Dimension] > node->GetUp()[*m_Dimension];
}



}// namespace itk

#endif
