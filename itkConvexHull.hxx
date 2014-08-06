#ifndef __itkConvexHull_hxx
#define __itkConvexHull_hxx

#include "itkConvexHull.h"

namespace itk
{

template< class TVector >
ConvexHull< TVector >
::ConvexHull()
{
	m_UpperHull = HalfConvexHullTreeContainerType::New();
	m_LowerHull = HalfConvexHullTreeContainerType::New();

	m_UpperHull->SetContour( Right );
	m_LowerHull->SetContour( Left );
}

template< class TVector >
bool ConvexHull< TVector >
::SetOtherDimension( const unsigned int dimension )
{
	if( dimension < Dimension )
	{
		m_OtherDimension = dimension;
		m_UpperHull->SetOtherDimension( m_OtherDimension );
		m_LowerHull->SetOtherDimension( m_OtherDimension );
		return true;
	}
	return false;
}

template< class TVector >
bool ConvexHull< TVector >
::SetDimension( const unsigned int dimension )
{
	if( dimension < Dimension )
	{
		m_Dimension = dimension;
		m_UpperHull->SetDimension( m_Dimension );
		m_LowerHull->SetDimension( m_Dimension );
		return true;
	}
	return false;
}

template< class TVector >
void ConvexHull< TVector >
::Add( VectorType & vector )
{
	if( m_UpperHull->Add( vector ) || m_LowerHull->Add( vector ) )
	{
		this->UpdateVerticalWidth( vector );
	}
}

template< class TVector >
void ConvexHull< TVector >
::Remove( VectorType & vector )
{
	if( m_UpperHull->Remove( vector ) || m_LowerHull->Remove( vector ) )
	{
		this->UpdateVerticalWidth();
	}
}

template< class TVector >
void ConvexHull< TVector >
::UpdateVerticalWidth()
{
	ConcatenableQueueNodeType *up = m_UpperHull->GetQRoot();
	ConcatenableQueueNodeType *low = m_LowerHull->GetQRoot(), tmp = low;
	double width = 0;
	while( up->HasChildren() )
	{
		ConcatenableQueueNodeType* node = up->GetChild( Right )->GetLowLeaf();
		if( tmp->GetDepth() < 3 )
		{
			m_VerticalWidth = 0;
			return;
		}
		while( tmp->GetDepth() > 2 )
		{
			if( *tmp < node )
			{
				tmp = tmp->GetChild( Right );
			}
			else
			{
				tmp = tmp->GetChild( Left );
			}
		}
		VectorType point = tmp->GetLow() + (double) ( node->Get()[m_Dimension] ) / ( tmp->GetUp()[m_Dimension] - tmp->GetLow()[m_Dimension] ) * ( tmp->GetUp() - tmp->GetLow() );
		double widthTmp = std::abs( node->Get()[m_OtherDimension] - point[m_OtherDimension] );
		if( widthTmp > width )
		{
			width = widthTmp;
			up = up->GetChild( Right );
		}
		else
		{
			up = up->GetChild( Left );
		}
	}
	m_VerticalWidth = width;
	width = 0;
	tmp = m_UpperHull->GetQRoot();
	while( low->HasChildren() )
	{
		ConcatenableQueueNodeType* node = low->GetChild( Right )->GetLowLeaf();
		if( tmp->GetDepth() < 3 )
		{
			m_VerticalWidth = 0;
			return;
		}
		while( tmp->GetDepth() > 2 )
		{
			if( *tmp < node )
			{
				tmp = tmp->GetChild( Right );
			}
			else
			{
				tmp = tmp->GetChild( Left );
			}
		}
		VectorType point = tmp->GetLow() + (double) ( node->Get()[m_Dimension] ) / ( tmp->GetUp()[m_Dimension] - tmp->GetLow()[m_Dimension] ) * ( tmp->GetUp() - tmp->GetLow() );
		double widthTmp = std::abs( node->Get()[m_OtherDimension] - point[m_OtherDimension] );
		if( widthTmp > width )
		{
			width = widthTmp;
			low = low->GetChild( Right );
		}
		else
		{
			low = low->GetChild( Left );
		}
	}
	if( width > m_VerticalWidth )
	{
		m_VerticalWidth = width;
	}
}

template< class TVector >
void ConvexHull< TVector >
::UpdateVerticalWidth( const VectorType & vector )
{
	ConcatenableQueueNodeType *right = m_UpperHull->GetQRoot(), *left = m_LowerHull->GetQRoot();
	ConcatenableQueueNodeType *node = right->GetNode( vector );
	ConcatenableQueueNodeType *notUpdated = left;
	if( node == NULL )
	{
		notUpdated = right;
	}
	if( notUpdated->GetDepth() < 3 )
	{
		m_VerticalWidth = 0;
		return;
	}
	while( notUpdated->GetDepth() > 2 )
	{
		if( *notUpdated < vector )
		{
			notUpdated = notUpdated->GetChild( Right );
		}
		else
		{
			notUpdated = notUpdated->GetChild( Left );
		}
	}
	VectorType point = notUpdated->GetLow() + (double) ( vector[m_Dimension] ) / ( notUpdated->GetUp()[m_Dimension] - notUpdated->GetLow()[m_Dimension] ) * ( notUpdated->GetUp() - notUpdated->GetLow() );
	double width = std::abs( vector[m_OtherDimension] - point[m_OtherDimension] );
	if( width > m_VerticalWidth )
	{
		m_VerticalWidth = width;
	}
}

}// namespace itk

#endif
