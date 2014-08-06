#ifndef __itkLineSegmentation_hxx
#define __itkLineSegmentation_hxx

#include "itkLineSegmentation.h"

#include <iostream>

namespace itk
{

template< class TIndexArray >
typename LineSegmentation< TIndexArray >
::DiscreteLineArrayPointer LineSegmentation< TIndexArray >
::SegmentLine()
{
	IndexConstIterator it = m_Input->begin();
	for( IndexConstIterator jt = it; jt != m_Input->end(); ++jt )
	{
		std::cout << *jt << " ";
	}
	std::cout << std::endl;
	this->Segment( m_Input->begin(), m_Input->end() );  
	this->SelectSegmentedLines();
	return m_SegmentedLine;
}

template< class TIndexArray >
void LineSegmentation< TIndexArray >
::SelectSegmentedLines()
{
	// TO DO 
}


template< class TIndexArray >
void LineSegmentation< TIndexArray >
::GCD( unsigned int a, unsigned int b ) const
{
	while( v != 0 )
	{
		unsigned int r = u % v;
		u = v;
		v = r;
	}
	return u;
}

template< class TIndexArray >
void LineSegmentation< TIndexArray >
::UpdateLine( ProjectedBlurredDiscreteLine2DType & line, const HalfConvexHullNodeType *LM, const HalfConvexHullNodeType *UM )
{
	unsigned int xdim = line.Low->GetDimension(), ydim = line.Low->GetOtherDimension();
	int r = line.MainVector[xdim] * LM->Get()[xdim] - line.MainVector[ydim] * LM->Get()[ydim];
	if( r == line.Offset )
	{
		line.U = UM;
	}
	else if( r == line.Offset + line.Thickness - 1 )
	{
		line.L = LM;
	}
	else if( r <= line.Offset - 1 )
	{
		HalfConvexHullNodeType *N = UM->GetLeft();
		UL = UM;
		int a0 = (M - N)[ydim], b0 = (M - N)[xdim];
		unsigned int gcd = this->GCD( a0, b0 );
		line.MainVector[xdim] = a0 / gcd;
		line.MainVector[ydim] = b0 / gcd;
		line.Offset = line.MainVector[xdim] * M[xdim] - line.MainVector[ydim] * M[ydim];
		while( ( line.L->Get()[xdim] - line.L->GetRight()->Get()[xdim] ) * line.MainVector[ydim] > 
			( line.L->Get()[ydim] - line.L->GetRight()->Get()[ydim] ) * line.MainVector[xdim] )
		{
			line.L = line.L->GetRight();
		}
		line.Thickness = line.MainVector[xdim] * line.L->Get()[xdim] - line.MainVector[ydim] * line.L->Get()[ydim] - line.Offset + 1;
	}
	else
	{
		HalfConvexHullNodeType *N = LM->GetLeft();
		line.L = LM;
		int a0 = (M - N)[ydim], b0 = (M - N)[xdim];
		unsigned int gcd = this->GCD( a0, b0 );
		line.MainVector[xdim] = a0 / gcd;
		line.MainVector[ydim] = b0 / gcd;
		line.Offset = line.MainVector[xdim] * M[xdim] - line.MainVector[ydim] * M[ydim];
		HalfConvexHullNodeType *node = line.L;
		while( ( line.U->Get()[xdim] - line.U->GetRight()->Get()[xdim] ) * line.MainVector[ydim] > 
			( line.U->Get()[ydim] - line.U->GetRight()->Get()[ydim] ) * line.MainVector[ydim] )
		{
			line.U = line.U->GetRight();
		}
		line.Thickness = line.MainVector[xdim] * line.L->Get()[xdim] - line.MainVector[ydim] * line.L->Get()[ydim] - line.Offset + 1;
	}
}

template< class TIndexArray >
float LineSegmentation< TIndexArray >
::GetBlurredLineWidth( const unsigned int omega, const int a, const int b ) const
{
	return ( float ) ( 1 - omega ) / ( std::abs( a ) > std::abs( b ) ? std::abs( a ) : std::abd( b ) );
}

template< class TIndexArray >
typename LineSegmentation< TIndexArray >
::LineType LineSegmentation< TIndexArray >
::InitializeLine( IndexConstIterator begin, unsigned int i, unsigned int j );
{
	LineType.Low = HalfConvexHullType::New();
	LineType.Up = HalfConvexHullType::New();
	LineType.Low->SetDimension( i );
	LineType.Low->SetOtherDimension( j );
	LineType.Up->SetDimension( i );
	LineType.Up->SetOtherDimension( j );
}

template< class TIndexArray >
::( IndexConstIterator begin, IndexConstIterator end )
{
	IndexArrayType MBS;
	IndexConstIterator C = begin;
	ProjectedDiscreteLine2DFixedArrayType lines;
	for( unsigned int i = 0, k = 0; i < Dimension; ++i )
	{
		for( unsigned int j = 0; j < i; ++j, ++k )
		{
			lines[k] = this->InitializeLine( C, i, j );
		}
	}
		
	unsigned int count = 0;
	do
	{
		count = 0;
		++C;
		for( unsigned int i = 0; i < Dimension; ++i )
		{
			lines[i].Low->Add( C );
			lines[i].Up->Add( C );
			if( this->UpdateLine( begin, C, lines[i] ) );
			{
				++count;
			}
		}
	} while( count > 1 );
	IndexConstIterator lastIndex = end;
	--lastIndex;
	MBS.push_back( IndexLine( begin,lastIndex ) );
	while( C != lastIndex )
	{
		do
		{
			count = 0;
			for( unsigned int i = 0; i < Dimension; ++i )
			{
				lines[i].Low->Remove( *begin );
				lines[i].Up->Remove( *begin );
			}
			++begin;
			for( unsigned int i = 0; i < Dimension; ++i )
			{
				if( !this->UpdateLine( begin, C, lines[i] ) );
				{
					++count;
				}
			}
		}while( count > 1 );

		do
		{
			count = 0;
			++C;
			for( unsigned int i = 0; i < Dimension; ++i )
			{
				lines[i].Low->Add( *C );
				lines[i].Up->Add( *C );
			}
			for( unsigned int i = 0; i < Dimension; ++i )
			{
				if( this->UpdateLine( begin, C, lines[i] ) );
				{
					++count;
				}
			}
		}while( count > 1 );
		lastIndex = C;
		--lastIndex;
		MBS.push_back( IndexLine( begin, lastIndex ) );
	}
	return MBS;
}

template< class TIndexArray >
void LineSegmentation< TIndexArray >
::UpdateLine3D( Projection2DBooleanFixedArrayType & lines, unsigned int id )
{
	if( id & 1 )
	{
		this->UpdateLine( lines[0], 
	}
}

template< class TIndexArray >
::( IndexConstIterator begin, IndexConstIterator end )
{
	IndexArrayType MBS;
	IndexConstIterator C = begin;
	unsigned int a = 0, b = 0, omega = b, mu = 0;
	bool x = this->GetBlurredLineWidth( e, a, b ) <= m_Width;
	bool y = this->GetBlurredLineWidth( e2, a, c ) <= m_Width;
	bool y = this->GetBlurredLineWidth( e, a, b ) <= m_Width;
	while( this->TwoBlurredSegmentInWidth( begin, C, m_Width ) )
	{
		++C;
		this->UpdateLine( begin, C, &a, &b, &c, &mu, &mu2, &e. &e2 );
	}
	IndexConstIterator lastIndex = end;
	--lastIndex;
	MBS.push_back( std::distance( begin, lastIndex ) );
	while( C != lastIndex )
	{
		while( this->TwoBlurredSegmentInWidth( begin, C, m_Width ) > m_Width )
		{
			++begin;
			this->UpdateLine( begin. C, &a, &b, &mu, &omega );
		}
		while( this->GetBlurredLineWidth( omega, a, b ) <= m_Width )
		{
			++C;
			this->UpdateLine( begin, C, &a, &b, &mu, &omega );
		}
		lastIndex = C;
		--lastIndex;
		MBS.push_back( std::distance( begin, lastIndex ) );
	}
}

templatE< class TIndexArray >
void LineSegmentation< TIndexArray >
::SegmentBlurredNaiveLine( IndexConstIterator begin, IndexConstIterator end ) const
{
	IndexConstIterator b = begin, e = begin;
	++e;
	Projection2DBooleanFixedArrayType lower;
	
	bool segment = true;
	unsigned int count = 0, minimumCount = 1;
	for( unsigned int i = 2; i < Dimension; ++i )
	{
		minimumCount *= i;
	}
	OffsetType offset = *begin - *e;
	ProjectedDiscreteLine2DFixedArrayType blurredNaiveLines2D;
	for( unsigned int i = 0, k = 0; i < Dimension; ++i )
	{
		for( unsigned int j = 0; j < Dimension; ++j, ++k )
		{
			blurredNaiveLines2D[k] = this->InitializeBlurredNaiveLine2D( b, e, i, j );
		}
	}
	++M;
	while( segment && M != end )
	{
		for( unsigned int i = 0; i < Projection2DNumber; ++i )
		{
			tmp[i] = this->ComputeExtendedBlurredNaiveLine2D( blurredNaiveLines2D[i], *e );
		}
		segment = this->IsBlurredNaiveLine( discreteLines2D );
		if( !segment )
		{
			for( unsigned int i = 0; i < Projection2DNumber; ++i )
			{
				discreteLines2D[i] = tmp[i];
			}
		}
		else
		{
			++e;
		}
	}
	MBS.push_back( this->ConvertTo3D( discreteLines2D, b, e ) )
	while( e != end )
	{
		do
		{
			++b;
			for( unsigned int i = 0; i < Projection2DNumber; ++i )
			{
				discreteLines2D[i] = this->ComputeCroppedBlurredNaiveLine( discreteLines2D[i], *b );
			}
		} while( !this->IsBlurredNaiveLine( discreteLines2D ) );
		do
		{
			++e;
			for( unsigned int i = 0; i < Projection2DNumber; ++i )
			{
				tmp[i] = this->ComputeExtendedBlurredNaiveLine2D( blurredNaiveLines2D[i], *e );
			}
		} while( this->IsBlurredNaiveLine( discreteLines2D ) );
		MBS.push_back( this->ConvertTo3D( discreteLines2D, b, e ) )
	}
	return MBS;
}

template< class TIndexArray >
LineSegmentation< TIndexArray >
::ProjectedBlurredDiscreteLine2DType LineSegmentation< TIndexArray >
::ComputeExtendedBlurredNaiveLine2D( const ProjectedBlurredDiscreteLine2DType & line, IndexConstIterator M, bool extension ) const
{
	ValueType r = this->R( line, *M );
	if( this->IsUpperLeaningPoint( line, r ) )
	{
		line.U = M;
	}
	if( this->IsLowerLeaningPoint( line, r ) )
	{
		line.L = M;
	}
	if( this->IsLowerExterior( line, r ) )
	{
		line.U = M;
		IndexConstIterator N = M;
		for( ; !this->IsUpperLeaningPoint( line, N ); --N );
		line.MainVector = this->Slope( N, M, line );
		line.Offset = this->R( line, *M );
		for( ; !this->IsLowerLeaning; );
	}
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsBlurredNaiveLine( const Projection2DFixedArrayType & discreteLines ) const
{
	unsigned int count;
	for( unsigned int i = 0; i < Projection2DNumber; ++i )
	{	
		if( this->IsBlurredNaiveLine2D( discreteLines2D[i] ) );
		{
			++count;
		}
	}
	return count < m_MinimumCount:
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsBlurredNaiveLine2D( const ProjectedDiscreteLine2DType & line ) const
{
	return line.Thickness - 1 <= std::trunc( line.Width * (double) line.MainVector[1] );
}

template< class TIndexArray >
double LineSegmentation< TIndexArray >
::VerticalWidth( const ProjectedBlurredDiscreteLine2DType & line ) const
{
	return (double) ( line.Thickness - 1 ) / (double) ( line.MainVector[1] );
}

template< class TIndexArray >
typename LineSegmentation< TIndexArray >
::ProjectedDiscreteLine2DType LineSegmentation< TIndexArray >
::InitializeDiscreteLine2D( IndexConstIterator begin, IndexConstIterator end, unsigned int x, unsigned int y )
{
	ProjectedDiscreteLine2DType line;
	line.Octant = this->GetOctant( *begin, *end );
	line.MainVector[0] = 0;
	line.MainVector[1] = 1;
	line,Offset = 0;
	line.Width = m_Width;
	line.X = x;
	line.Y = y;
	line.Thickness = 1;
	line.L = *begin;
	line.U = *begin;
	line.Begin = begin;
	line.End = end;
}

template< class TIndexArray >
void LineSegmentation< TIndexArray >
::Segment( IndexConstIterator begin, IndexConstIterator end )
{
	m_SegmentedLine->clear();
	IndexConstIterator M = begin;
	IndexConstIteratorContainerType UF, LF, UL, LL;
	++M;
	for( unsigned int i = 0; i < Dimension; ++i )
	{
		UF.push_back( begin );
		LF.push_back( begin );
		UL.push_back( M );
		LL.push_back( M );
	}
		
}

template< TIndexArray >
void LineSegmentation< TIndexArray >
::DiscreteLine2DOctantTransform( DiscreteLine2DType & line, unsigned int octant ) const
{
	line.MainVector = this->OctantTransform< Vector2DType >( line.MainVector, octant );
	if( octant > 0 && octant < 5 )
	{
		line.Offset = -line.Offset - line.Thickness + 1;
	}
	if( octant > 2 && octant < 7 )
	{
		ValueType tmp = -transform[0];
		transform[0] = -transform[1];
		transform[1] = tmp;	
	}
}

template< TIndexArray >
void LineSegmentation< TIndexArray >
::FirstBissectriceSymmetry( DiscreteLine2DType & line, bool increasing ) const
{
	if( !increasing )
	{
		ValueType tmp = -line.MainVector[0];
		line.MainVector[0] = -line.MainVector[1];
		line.MainVector[1] = tmp;
	}
	else
	{
		ValueType tmp = line.MainVector[0];
		line.MainVector[0] = line.MainVector[1];
		line.MainVector[1] = tmp;
		line.Offset = -line.Offset - line.Thickness + 1;
	}
}

template< TIndexArray >
void LineSegmentation< TIndexArray >
::OxSymmetry( DiscreteLine2DType & line, bool increasing ) const
{
	if( !increasing )
	{
		line.MainVector[1] = -line.MainVector[1];
	}
	else
	{
		line.MainVector[0] = -line.MainVector[0];
		line.Offset = -line.Offset - line.Thickness + 1;
	}
}

template< TIndexArray >
void LineSegmentation< TIndexArray >
::OySymmetry( DiscreteLine2DType & line, bool increasing ) const
{
	if( increasing )
	{
		line.MainVector[0] = -line.MainVector[0];
	}
	else
	{
		line.MainVector[1] = -line.MainVector[1];
		line.Offset = -line.Offset - line.Thickness + 1;
	}
}

template< TIndexArray >
void LineSegmentation< TIndexArray >
::OSymmetry( DiscreteLine2DType & line, bool increasing ) const
{
	if( !increasing )
	{
		line.MainVector[0] = -line.MainVector[0];
		line.MainVector[1] = -line.MainVector[1];
	}
	else
	{
		line.Offset = -line.Offset - line.Thickness + 1;
	}
}

template< TIndexArray >
void LineSegmentation< TIndexArray >
::OctantTransform( DiscreteLine2DType & line, unsigned int octant, bool increasing ) const
{
	switch( octant )
	{
		case 0:
			break;
		case 1:
			this->FirstBissectriceSymmetry( line, increasing );
			break;
		case 2:
			this->OySymmetry( line );
			this->FirstBissectriceSymmetry( line, increasing );
			break;
		case 3:
			this->OxSymmetry( line, increasing );
			break;
		case 4:
			this->OSymmetry( line, increasing );
			break;
		case 5:
			this->OSymmetry( line, increasing );
			this->FirstBissectriceSymmetry( line, increasing );
			break;
		case 6:
			this->OxSymmetry( line, increasing );
			this->FirstBissectriceSymmetry( line, increasing );
			break;
		case 7:
			this->OxSymmetry( line, increasing );

	}
}

template< TIndexArray >
template< TVector >
unsigned int LineSegmentation< TIndexArray >
::GetOctant( const TVector & vector, unsigned int i, unsigned int j )
{
	int x = vector[i], y = vector[j];
	unsigned int octant = 0;
	if( y > 0 )
	{   
		octant |= 4;
	}   
	if( y * x < 0 )
	{
		octant |= 2;
		if( std::abs( y ) > std::abs( x ) )
		{   
			octant |= 1;
		}   
	}
	else
	{   
		if( std::abs( y ) < std::abs( x ) )
		{   
			octant |= 1;
		}   
	}
	return octant;
}   

template< class TIndexArray >
void LineSegmentation< TIndexArray >
::Segment( IndexConstIterator begin, IndexConstIterator end )
{
	m_SegmentedLine->clear();
	bool lower = true;
	IndexConstIterator M = begin;
	IndexConstIteratorContainerType UF, LF, UL, LL;
	++M;
	for( unsigned int i = 0; i < Dimension; ++i )
	{
		UF.push_back( begin );
		LF.push_back( begin );
		UL.push_back( M );
		LL.push_back( M );
	}
	unsigned int dimension = 0;
	DiscreteLineType line = this->ComputeDiscreteLine( LF, UF, LL, UL, begin, M, lower, dimension );
	std::cout << " line " << " MainVector " << line.MainVector << " Offset " << line.Offset << std::endl;
	while( M != end )
	{
		bool keepGoing = true;
		while( M != end && keepGoing )
		{
			std::cout << " M " << *M << std::endl;
			DiscreteLineType tmp = line;
			while( this->IsWeaklyExterior( *M, line, lower, dimension ) )
			{
				std::cout << " Is Weakly Exterior " << std::endl;
				tmp = this->ComputeDiscreteLine( LF, UF, LL, UL, begin, M, lower, dimension );
				std::cout << " line " << " MainVector " << line.MainVector << " Offset " << line.Offset << std::endl;
				this->UpdateLU( LF, UF, LL, UL, line, begin, M );
			}
			while( this->IsStronglyExterior( *M, line, lower, dimension ) )
			{
				std::cout << " Is Strongly Exterior " << std::endl;
				m_SegmentedLine->push_back( line );
				keepGoing = false;
				do
				{
					std::cout << " begin " << *begin << std::endl;
					++begin;
				} while( this->IsExterior( *begin, tmp ) );
				tmp = this->ComputeDiscreteLine( LF, UF, LL, UL, begin, M, lower, dimension );
				std::cout << " line " << " MainVector " << line.MainVector << " Offset " << line.Offset << std::endl;
				this->UpdateLU( LF, UF, LL, UL, line, begin, M );
			}
			line = tmp;
			++M;
		}
	}		

	DiscreteLineType tmp = line;
	for( unsigned int i = 0; i < Dimension; ++i )
	{
		if( i != line.MainDimension )
		{
			lower = this->IsWeaklyExteriorLower( *M, tmp, i );
			if( lower || this->IsWeaklyExteriorUpper( *M, tmp, i ) )
			{
				this->ComputeDiscreteLine( tmp, LF, UF, LL, UL, begin, M, lower, i );
				this->UpdateLU( LF, UF, LL, UL, tmp, begin, M, i );
			}
		}
	}
	for( unsigned int i = 0; i < Dimension; ++i )
	{
		if( i != line.MainDimension )
		{
			lower = this->IsStronglyExteriorLower( *M, tmp, i );
			if( lower || this->IsStronglyExteriorUpper( *M, tmp, i ) )
			{
				m_SegmentedLine->push_back( line );
				bool keepGoing = false;
				do
				{
					++begin;
					for( unsigned int j = 0; j < Dimension; ++j )
					{
						if( j != line.MainDimension && this->IsExterior( *begin, tmp, j ) )
						{
							keepGoing = true;
						}
					}
				}while( keepGoing );
				
			}
		}
	}
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::UpdateLU( IndexConstIteratorContainerType & LF, IndexConstIteratorContainerType & UF, IndexConstIteratorContainerType & LL, IndexConstIteratorContainerType & UL, const DiscreteLineType & line, IndexConstIterator it, IndexConstIterator jt ) const
{
	std::vector< bool > LProcessed, UProcessed;
	for( unsigned int i  = 0; i < Dimension; ++i )
	{
		LProcessed.push_back( false );
		UProcessed.push_back( false );
	}
	IndexConstIterator kt = jt;
	do
	{
		for( unsigned int i = 0; i < Dimension; ++i )
		{
			if( i != line.MainDimension )
			{
				if( this->IsOnLowerBoundary( *kt, line, i ) )
				{
					if( !LProcessed[i] )
					{
						LL[i] = kt;
						LProcessed[i] = true;
					}
					LF[i] = kt;
				}
				if( this->IsOnUpperBoundary( *kt, line, i ) )
				{
					if( !UProcessed[i] )
					{
						UL[i] = kt;
						UProcessed[i] = true;
					}
					UF[i] = kt;
				}
			}
		}
	} while( kt-- != it );
	return true;
}

template< class TIndexArray >
typename LineSegmentation< TIndexArray >
::DiscreteLineType LineSegmentation< TIndexArray >
::ComputeDiscreteLine( IndexConstIteratorContainerType & LF, IndexConstIteratorContainerType & UF, IndexConstIteratorContainerType & LL, IndexConstIteratorContainerType & UL, IndexConstIterator begin, IndexConstIterator it, const std::vector< bool > & lower, unsigned int & dimension ) const
{
	IndexType XF;
	if( lower[dimension] )
	{
		std::cout << " lower " << std::endl;
		XF = *( LF[dimension] );
	}
	else
	{
		std::cout << " upper " << std::endl;
		XF = *( UF[dimension] );
	}
	std::cout << " Compute Discrete Line between " << *it << " and " << *begin << std::endl;
	std::cout << " XF " << XF << std::endl;
	OffsetType vector = *it - XF;
	DiscreteLineType line;
	ValueType max = 0;
	unsigned int imax = 0;
	for( unsigned int i = 0; i < Dimension; ++i )
	{
		line.MainVector[i] =  vector[i];
		if( std::abs( vector[i] ) > max )
		{
			max = std::abs( vector[i] );
			imax = i;
		}
	}
	line.MainDimension = imax;
	for( unsigned int i = 0; i < Dimension; ++i )
	{
		line.Thickness[i] = std::abs( line.MainVector[line.MainDimension] );
	}
	for( unsigned int i = 0; i < Dimension; ++i )
	{
		if( i == line.MainDimension )
		{
			line.Offset[i] =  0;
		}
		else
		{
			if( lower[i] )
			{	
				line.Offset[i] = this->Cost( *( UL[dimension] ), line, dimension );
			}
			else
			{
				line.Offset[i] = this->Cost( *( LL[dimension] ), line, dimension ) - line.Thickness[i] + 1;
			}
		}
	}
	line.Begin = begin;
	return line;
}

template< class TIndexArray >
typename LineSegmentation< TIndexArray >
::ValueType LineSegmentation< TIndexArray >
::Cost( const Index2DType & idx, const Vector2DType & vector ) const
{
	return vector[1] * idx[0] - vector[0] * idx[1];
}

template< class TIndexArray >
typename LineSegmentation< TIndexArray >
::ValueType LineSegmentation< TIndexArray >
::Cost( const Index2DType & idx, const DiscreteLine2DType & line ) const
{
	return line.MainVector[1] * idx[0] - line.MainVector[0] * idx[1];
}

template< class TIndexArray >
typename LineSegmentation< TIndexArray >
::ValueType LineSegmentation< TIndexArray >
::Cost( const IndexType & idx, const DiscreteLineType & line, unsigned int i ) const
{
	return line.MainVector[i] * idx[line.MainDimension] - line.MainVector[line.MainDimension] * idx[i];
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsIncreasing( const DiscreteLineType & line, unsigned int & i ) const
{
	return line.MainVector[line.MainDimension] * line.MainVector[i] >= 0;
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsDecreasing( const DiscreteLineType & line, unsigned int & i ) const
{
	return line.MainVector[line.MainDimension] * line.MainVector[i] <= 0;
}



template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsOnBoundary( const IndexType & M, const DiscreteLineType & line ) const
{
	for( unsigned int  i = 0; i < Dimension; ++i )
	{
		if( i != line.MainDimension )
		{
			bool tmp = this->IsOnLowerBoundary( M, line, i );
			if( tmp || this->IsOnUpperBoundary( M, line, i ) ) 
			{
				return true;
			}
		}
	}
	return false;;
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsOnBoundary( const IndexType & M, const DiscreteLineType & line, bool & lower, unsigned int & dimension ) const
{
	for( unsigned int  i = 0; i < Dimension; ++i )
	{
		if( i != line.MainDimension )
		{
			bool tmp = this->IsOnLowerBoundary( M, line, i );
			if( tmp || this->IsOnUpperBoundary( M, line, i ) ) 
			{
				dimension = i;
				lower = tmp;
				return true;
			}
		}
	}
	return false;
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsExterior( const IndexType & M, const DiscreteLineType & line ) const
{
	for( unsigned int  i = 0; i < Dimension; ++i )
	{
		if( i != line.MainDimension )
		{
			bool tmp = this->IsExteriorLower( M, line, i );
			if( tmp || this->IsExteriorUpper( M, line, i ) ) 
			{
				return true;
			}
		}
	}
	return false;
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsExterior( const IndexType & M, const DiscreteLineType & line, bool & lower, unsigned int & dimension ) const
{
	for( unsigned int  i = 0; i < Dimension; ++i )
	{
		if( i != line.MainDimension )
		{
			bool tmp = this->IsExteriorLower( M, line, i );
			if( tmp || this->IsExteriorUpper( M, line, i ) ) 
			{
				dimension = i;
				lower = tmp;
				return true;
			}
		}
	}
	return false;
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsWeaklyExterior( const IndexType & M, const DiscreteLineType & line ) const
{
	for( unsigned int  i = 0; i < Dimension; ++i )
	{
		if( i != line.MainDimension )
		{
			bool tmp = this->IsWeaklyExteriorLower( M, line, i );
			if( tmp || this->IsWeaklyExteriorUpper( M, line, i ) ) 
			{
				return true;
			}
		}
	}
	return false;
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsWeaklyExterior( const IndexType & M, const DiscreteLineType & line, bool & lower, unsigned int & dimension ) const
{
	for( unsigned int  i = 0; i < Dimension; ++i )
	{
		if( i != line.MainDimension )
		{
			bool tmp = this->IsWeaklyExteriorLower( M, line, i );
			if( tmp || this->IsWeaklyExteriorUpper( M, line, i ) ) 
			{
				dimension = i;
				lower = tmp;
				return true;
			}
		}
	}
	return false;
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsStronglyExterior( const IndexType & M, const DiscreteLineType & line ) const
{
	for( unsigned int  i = 0; i < Dimension; ++i )
	{
		if( i != line.MainDimension )
		{
			bool tmp = this->IsStronglyExteriorLower( M, line, i );
			if( tmp || this->IsStronglyExteriorUpper( M, line, i ) ) 
			{
				return true;
			}
		}
	}
	return false;
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsStronglyExterior( const IndexType & M, const DiscreteLineType & line, bool & lower, unsigned int & dimension ) const
{
	for( unsigned int  i = 0; i < Dimension; ++i )
	{
		if( i != line.MainDimension )
		{
			bool tmp = this->IsStronglyExteriorLower( M, line, i );
			if( tmp || this->IsStronglyExteriorUpper( M, line, i ) ) 
			{
				dimension = i;
				lower = tmp;
				return true;
			}
		}
	}
	return false;
}


template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsExteriorLower( const ValueType & value, const DiscreteLine2DType & line ) const
{
	return value >= ( line.Offset + line.Thickness );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsExteriorLower( const Index2DType & idx, const DiscreteLine2DType & line ) const 
{
	return this->IsExteriorLower( this->Cost( idx, line ), line );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsExteriorUpper( const ValueType & value, const DiscreteLine2DType & line ) const 
{
	return value <= line.Offset - 1;
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsExteriorUpper( const Index2DType & idx, const DiscreteLine2DType & line ) const 
{
	return this->IsWeaklyExteriorUpper( this->Cost( idx, line ), line );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsExterior( const ValueType & value, const DiscreteLine2DType & line ) const 
{
	return this->IsExteriorUpper( value, line ) || this->IsExteriorLower( value, line );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsExterior( const Index2DType & idx, const DiscreteLine2DType & line ) const 
{
	ValueType value = this->Cost( idx, line );
	return this->IsExterior( value, line );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsWeaklyExteriorLower( const ValueType & value, const DiscreteLine2DType & line ) const
{
	return value == ( line.Offset + line.Thickness );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsWeaklyExteriorLower( const Index2DType & idx, const DiscreteLine2DType & line ) const 
{
	return this->IsWeaklyExteriorLower( this->Cost( idx, line ), line );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsWeaklyExteriorUpper( const ValueType & value, const DiscreteLine2DType & line ) const 
{
	return value == line.Offset - 1;
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsWeaklyExteriorUpper( const Index2DType & idx, const DiscreteLine2DType & line ) const 
{
	return this->IsWeaklyExteriorUpper( this->Cost( idx, line ), line );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsWeaklyExterior( const ValueType & value, const DiscreteLine2DType & line ) const 
{
	return this->IsWeaklyExteriorUpper( value, line ) || this->IsWeaklyExteriorLower( value, line );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsOnLowerBoundary( const ValueType & value, const DiscreteLine2DType & line ) const
{
	return value == ( line.Offset + line.Thickness - 1 );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsOnLowerBoundary( const Index2DType & idx, const DiscreteLine2DType & line ) const 
{
	return this->IsOnLowerBoundary( this->Cost( idx, line ), line );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsOnUpperBoundary( const ValueType & value, const DiscreteLine2DType & line ) const 
{
	return value == line.Offset;
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsOnUpperBoundary( const Index2DType & idx, const DiscreteLine2DType & line ) const 
{
	return this->IsOnUpperBoundary( this->Cost( idx, line ), line );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsOnBoundary( const ValueType & value, const DiscreteLine2DType & line ) const 
{
	return this->IsOnUpperBoundary( value, line ) || this->IsOnLowerBoundary( value, line );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsWeaklyExterior( const Index2DType & idx, const DiscreteLine2DType & line ) const 
{
	ValueType value = this->Cost( idx, line );
	return this->IsWeaklyExterior( value, line );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsStronglyExteriorLower( const ValueType & value, const DiscreteLine2DType & line ) const
{
	return value > ( line.Offset + line.Thickness );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsStronglyExteriorLower( const Index2DType & idx, const DiscreteLine2DType & line ) const 
{
	return this->IsStronglyExteriorLower( this->Cost( idx, line ), line );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsStronglyExteriorUpper( const ValueType & value, const DiscreteLine2DType & line ) const 
{
	return value < line.Offset - 1;
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsStronglyExteriorUpper( const Index2DType & idx, const DiscreteLine2DType & line ) const 
{
	return this->IsStronglyExteriorUpper( this->Cost( idx, line ), line );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsStronglyExterior( const ValueType & value, const DiscreteLine2DType & line ) const 
{
	return this->IsStronglyExteriorUpper( value, line ) || this->IsStronglyExteriorLower( value, line );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsStronglyExterior( const Index2DType & idx, const DiscreteLine2DType & line ) const 
{
	ValueType value = this->Cost( idx, line );
	return this->IsStronglyExterior( value, line );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsExteriorLower( const ValueType & value, const DiscreteLineType & line, unsigned int i ) const
{
	return value >= ( line.Offset[i] + line.Thickness[i] );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsExteriorLower( const IndexType & idx, const DiscreteLineType & line, unsigned int i ) const 
{
	return this->IsExteriorLower( this->Cost( idx, line, i ), line, i );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsExteriorUpper( const ValueType & value, const DiscreteLineType & line, unsigned int i ) const 
{
	return value <= line.Offset[i] - 1;
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsExteriorUpper( const IndexType & idx, const DiscreteLineType & line, unsigned int i ) const 
{
	return this->IsExteriorUpper( this->Cost( idx, line, i ), line, i );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsExterior( const ValueType & value, const DiscreteLineType & line, unsigned int i ) const 
{
	return this->IsExteriorUpper( value, line, i ) || this->IsExteriorLower( value, line, i );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsExterior( const IndexType & idx, const DiscreteLineType & line, unsigned int i ) const 
{
	ValueType value = this->Cost( idx, line, i );
	return this->IsExterior( value, line, i );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsWeaklyExteriorLower( const ValueType & value, const DiscreteLineType & line, unsigned int i ) const
{
	return value == ( line.Offset[i] + line.Thickness[i] );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsWeaklyExteriorLower( const IndexType & idx, const DiscreteLineType & line, unsigned int i ) const 
{
	return this->IsWeaklyExteriorLower( this->Cost( idx, line, i ), line, i );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsWeaklyExteriorUpper( const ValueType & value, const DiscreteLineType & line, unsigned int i ) const 
{
	return value == line.Offset[i] - 1;
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsWeaklyExteriorUpper( const IndexType & idx, const DiscreteLineType & line, unsigned int i ) const 
{
	return this->IsWeaklyExteriorUpper( this->Cost( idx, line, i ), line, i );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsWeaklyExterior( const ValueType & value, const DiscreteLineType & line, unsigned int i ) const 
{
	return this->IsWeaklyExteriorUpper( value, line, i ) || this->IsWeaklyExteriorLower( value, line, i );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsWeaklyExterior( const IndexType & idx, const DiscreteLineType & line, unsigned int i ) const 
{
	ValueType value = this->Cost( idx, line, i );
	return this->IsWeaklyExterior( value, line, i );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsOnLowerBoundary( const ValueType & value, const DiscreteLineType & line, unsigned int i ) const
{
	return value == ( line.Offset[i] + line.Thickness[i] - 1);
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsOnLowerBoundary( const IndexType & idx, const DiscreteLineType & line, unsigned int i ) const 
{
	return this->IsOnLowerBoundary( this->Cost( idx, line, i ), line, i );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsOnUpperBoundary( const ValueType & value, const DiscreteLineType & line, unsigned int i ) const 
{
	return value == line.Offset[i];
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsOnUpperBoundary( const IndexType & idx, const DiscreteLineType & line, unsigned int i ) const 
{
	return this->IsOnUpperBoundary( this->Cost( idx, line, i ), line, i );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsOnBoundary( const ValueType & value, const DiscreteLineType & line, unsigned int i ) const 
{
	return this->IsOnUpperBoundary( value, line, i ) || this->IsOnLowerBoundary( value, line, i );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsOnBoundary( const IndexType & idx, const DiscreteLineType & line, unsigned int i ) const 
{
	ValueType value = this->Cost( idx, line, i );
	return this->IsOnBoundary( value, line, i );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsStronglyExteriorLower( const ValueType & value, const DiscreteLineType & line, unsigned int i ) const
{
	return value > ( line.Offset[i] + line.Thickness[i] );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsStronglyExteriorLower( const IndexType & idx, const DiscreteLineType & line, unsigned int i ) const 
{
	return this->IsStronglyExteriorLower( this->Cost( idx, line, i ), line, i );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsStronglyExteriorUpper( const ValueType & value, const DiscreteLineType & line, unsigned int i ) const 
{
	return value < line.Offset[i] - 1;
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsStronglyExteriorUpper( const IndexType & idx, const DiscreteLineType & line, unsigned int i ) const 
{
	return this->IsStronglyExteriorUpper( this->Cost( idx, line, i ), line, i );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsStronglyExterior( const ValueType & value, const DiscreteLineType & line, unsigned int i ) const 
{
	return this->IsStronglyExteriorUpper( value, line, i ) || this->IsStronglyExteriorLower( value, line, i );
}

template< class TIndexArray >
bool LineSegmentation< TIndexArray >
::IsStronglyExterior( const IndexType & idx, const DiscreteLineType & line, unsigned int i ) const 
{
	ValueType value = this->Cost( idx, line );
	return this->IsStronglyExterior( value, line );
}

}// namespace itk

#endif
