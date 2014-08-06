#ifndef __itkBlurred3DLineSegmentation_hxx
#define __itkBlurred3DLineSegmentation_hxx

#include "itkBlurred3DLineSegmentation.h"

namespace itk
{

template< class TIndexArray >
template< class TIterator >
typename Blurred3DLineSegmentation< TIndexArray >
::MaxDiscreteSegmentContainerType Blurred3DLineSegmentation< TIndexArray >
::Segment( TIterator begin, TIterator end )
{
	DiscreteLineType line;
	MaxDiscreteSegmentContainerType MBS;
//	DiscreteSegmentType segment;
	typename ConvexHullType::Pointer hull[3];
	for( unsigned int i = 0; i < Dimension; ++i )
	{
		hull[i] = ConvexHullType::New();
	}
	this->Add( hull, *begin );
	unsigned int projId = 0, count = 0;
	TIterator current = begin;
	++current;
	do
	{
		count = 0;
		for( unsigned int i = 0; i < Dimension; ++i )
		{
			if( this->IsBlurredLine( line[i] ) )
			{
				projId |= ( 1 << i );
				++count;
			}
		}
		if( count > 1 )
		{
			this->Add( hull, *current );
			this->UpdateLine( line, hull );
			++current;
		}
	}while( count > 1 );

//	segment.Begin = begin;
//	segment.End-- = current;
//	segment.Line = line;
//	segment.Dim = projId;
//	MBS.push_back( segment );

	while( current != end )
	{
		do
		{
			count = 0;
			for( unsigned int i = 0; i < Dimension; ++i )
			{
				if( !this->IsBlurredLine( line[i] ) )
				{
					projId |= ( 1 << i );
					++count;
				}
			}
			if( count > 1 )
			{
				this->Remove( hull, *begin );
				this->UpdateLine( line, hull );
				++begin;
			}
		}while( count > 1 && current != end );
		do
		{
			count = 0;
			for( unsigned int i = 0; i < Dimension; ++i )
			{
				if( !this->IsBlurredLine( line[i] ) )
				{
					projId |= ( 1 << i );
					++count;
				}
			}
			if( count > 1 )
			{
				this->Add( hull, *current );
				this->UpdateLine( line, hull );
				++current;
			}
		}while( count > 1 && current != end );
//		segment.Begin = begin;
//		segment.End-- = current;
//		segment.Line = line;
//		segment.Dim = projId;
	//	MBS.push_back( segment );
	}
	return MBS;
}

template< class TIndexArray >
void Blurred3DLineSegmentation< TIndexArray >
::UpdateLine( const ConvexHullType *hull, Discrete2DLineType & line, ConcatenableQueueNodeType *M, bool increase ) const
{
	if( increase )
	{
		int r = this->GetRemainder( M, line );
		if( r == line.Offset )
		{
			line.UL = M;
		}
		else if( r == line.Offset + line.Width - 1 )
		{
			line.LL = M;
		}
		else if( r <= line.Offset - 1 )
		{
			line.UL = M;
			VectorType N = *(hull->GetUpperHull()->GetQRoot()->GetUpLeaf()->GetLeft());
			int a = M->Get()[line.OtherDimension] - N[line.OtherDimension];
			int b = M->Get()[line.Dimension] - N[line.Dimension];
			ConcatenableQueueNodeType *C = line.LL;
			while( this->Slope( *C, *(C->GetRight()) ) <= (double) (line.a) / (double) (line.b) )
			{
				C = C->GetRight();
			}
			line.LL = C;
		}
		else
		{
			line.LL = M;
			VectorType N = *(hull->GetLowerHull()->GetQRoot()->GetUpLeaf()->GetLeft());
			int a = M->Get()[line.OtherDimension] - N[line.OtherDimension];
			int b = M->Get()[line.Dimension] - N[line.Dimension];
			ConcatenableQueueNodeType *C = line.LL;
			while( this->Slope( *C, *(C->GetRight()) ) <= (double) (line.a) / (double) (line.b) )
			{
				C = C->GetRight();
			}
			line.LL = C;
		}
	}
	else
	{
		// TODO decrease case
	}
	line.Width = hull->GetVerticalWidth();

}

template< class TIndexArray >
bool Blurred3DLineSegmentation< TIndexArray >
::IsBlurredLine( const Discrete2DLineType & line ) const
{
	return line.Thickness - 1 < m_Width * line.b;
}

template< class TIndexArray >
void Blurred3DLineSegmentation< TIndexArray >
::Add( ConvexHullType *hull, VectorType point ) const
{
	hull->Add( point );
}

template< class TIndexArray >
void Blurred3DLineSegmentation< TIndexArray >
::Remove( ConvexHullType *hull, VectorType point ) const
{
	hull->Remove( point );
}

} // namespace itk

#endif
