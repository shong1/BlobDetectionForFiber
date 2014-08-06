#ifndef __itkSkeletonToGraphFilter_hxx
#define __itkSkeletonToGraphFilter_hxx

#include "itkSkeletonToGraphFilter.h"

#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodIterator.h"

#include <iostream>

namespace itk
{

// Constructor
template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::SkeletonToGraphFilter()
{
	this->ProcessObject::SetNumberOfRequiredInputs( 2 );

	m_NeighborhoodCount = 1;

	for( unsigned int i = 0; i < ImageDimension; ++i )
	{
		m_NeighborhoodCount *= 3;
	}

	m_Label = 1;

	GraphPointer output = dynamic_cast< GraphType* > ( this->MakeOutput( 0  ).GetPointer() );
	LabelImagePointer outputLabel = dynamic_cast< LabelImageType* > ( this->MakeOutput( 1  ).GetPointer() );
	IdentifierImagePointer outputIdentifier = dynamic_cast< IdentifierImageType* > ( this->MakeOutput( 2  ).GetPointer() );

	this->ProcessObject::SetNumberOfRequiredOutputs( 3 );
	this->ProcessObject::SetNthOutput( 1, outputLabel.GetPointer() );
	this->ProcessObject::SetNthOutput( 2, outputIdentifier.GetPointer() );
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::~SkeletonToGraphFilter() 
{}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
DataObject::Pointer
SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::MakeOutput( unsigned int idx )
{
	switch( idx )
	{
		case 0: {
			GraphPointer outputGraph = GraphType::New();
			return dynamic_cast< DataObject * > ( outputGraph.GetPointer() );
			break; }
		case 1: {
			LabelImagePointer outputLabel = LabelImageType::New();
			return dynamic_cast< DataObject * > ( outputLabel.GetPointer() );
			break; }
		case 2: {
			IdentifierImagePointer outputIdentifier = IdentifierImageType::New();
			return dynamic_cast< DataObject * > ( outputIdentifier.GetPointer() );
			break; }
		default:
			return NULL;
			break;
	}
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
void SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::SetInput( ImageType * image )
{
	this->SetNthInput( 0, const_cast< ImageType * > ( image ) );
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
void SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::SetInput( unsigned int idx, const ImageType * input )
{
	this->ProcessObject::SetNthInput( idx, const_cast< ImageType * > ( input ) );
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
const typename SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::ImageType * SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::GetInput( unsigned int idx )
{
	return dynamic_cast< const ImageType * > ( this->ProcessObject::GetInput( idx ) );
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
typename SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::GraphType * SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::GetOutput( void )
{
	return dynamic_cast< GraphType * > ( this->ProcessObject::GetOutput( 0 ) );
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
typename SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::LabelImageType * SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::GetLabelOutput( void )
{
	return dynamic_cast< LabelImageType * > ( this->ProcessObject::GetOutput( 1 ) );
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
typename SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::IdentifierImageType * SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::GetIdentifierOutput( void )
{
	return dynamic_cast< IdentifierImageType * > ( this->ProcessObject::GetOutput( 2 ) );
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
void SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::GenerateData()
{
	ImageConstPointer input = static_cast< const ImageType * > ( this->GetInput( 0 ) );

	typename ImageIteratorType::RadiusType radius;
	radius.Fill( 1 );

	// Initialize every images

	this->GetLabelOutput()->SetRegions( input->GetLargestPossibleRegion() );
	this->GetLabelOutput()->Allocate();
	this->GetLabelOutput()->FillBuffer( 0 );
	this->GetLabelOutput()->Update();

	this->GetIdentifierOutput()->SetRegions( input->GetLargestPossibleRegion() );
	this->GetIdentifierOutput()->Allocate();
	this->GetIdentifierOutput()->FillBuffer( 0 );
	this->GetIdentifierOutput()->Update();

	this->GetOutput()->GetEdgeContainer()->Reserve( 1000000 );

	// Initialize iterators

	ImageIteratorType it( radius, input, input->GetLargestPossibleRegion() );
	LabelImageIteratorType lit( radius, this->GetLabelOutput(), input->GetLargestPossibleRegion() );
	IdentifierImageIteratorType iit( radius, this->GetIdentifierOutput(), input->GetLargestPossibleRegion() );
	VesselnessImageIteratorType vit( radius, this->GetVesselnessInput(), input->GetLargestPossibleRegion() );

	m_Region = input->GetLargestPossibleRegion();

	GraphPointer output = this->GetOutput();

	LabelType label = 1;

	unsigned int neighborhoodCount = std::pow( 3, ImageDimension );
	unsigned int center = neighborhoodCount / 2;

	it.GoToBegin();
	lit.GoToBegin();
	iit.GoToBegin();
	vit.GoToBegin();

	// Go through the whole image
	while( this->ReachBranchSkeleton( it, lit, iit, vit ) )
	{
		ImageIteratorType tit = it;
		LabelImageIteratorType tlit = lit;
		IdentifierImageIteratorType tiit = iit;
		VesselnessImageIteratorType tvit = vit;

		this->TransformSkeletonToGraph( tit, tlit, tiit, tvit );
	}
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
void SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::TransformSkeletonToGraph( ImageIteratorType & it, LabelImageIteratorType & lit, IdentifierImageIteratorType & iit, VesselnessImageIteratorType & vit )
{
	NodePointerType node = this->GetOutput()->CreateNewNode();
	node->LineIndexes = IndexContainerType();

	node->Weight = 1;
	node->BackLabel = 0;
	node->FrontLabel = 0;

	NodePointerType nodes[2];

	this->ConnectSelf( node );

	ImageIteratorType tit = it;
	LabelImageIteratorType tlit = lit;
	IdentifierImageIteratorType tiit = iit;
	VesselnessImageIteratorType tvit = vit;


	bool borderFound = false;
	this->MakeNode( tit, tlit, tiit, tvit, node, &borderFound );
	// If the line is not fully scanned, we recall the function
	if( !borderFound )
	{
		borderFound = true;
		this->MakeNode( it, lit, iit, vit, node, &borderFound );
	}
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
void SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::MakeNode( ImageIteratorType & it, LabelImageIteratorType & lit, IdentifierImageIteratorType & iit, VesselnessImageIteratorType & vit, NodePointerType node, bool * borderFound )
{ 
	ArrayType map;
	unsigned int count = this->GetNeighborhoodMap( it, map );

	bool end = false;

	// If we are not in the middle of a line
	if( count != 2 )
	{

		// It is an intersection
		if( count > 2 )
		{
			ImageIteratorType tit = it;
			LabelImageIteratorType tlit = lit;
			IdentifierImageIteratorType tiit = iit;
			VesselnessImageIteratorType tvit = vit;
			
			// Node not existing -> creating new node
			// The label image helps to know which line has been scanned
			if( !lit.GetCenterPixel() )
			{
				this->FillIntersection( tit, tlit, tiit, tvit, m_Label );
			}
			// Already scanned intersection
			else
			{
				this->ConnectNodes( tit, tlit, tiit, tvit, node, m_Label );
				// We increment label to not separate every intersection
				++m_Label;
			}
		}

		// If we already found a border, and we are not anymore in the
		// middle of the line, this is the end
		if( *borderFound )
		{
			end = true;
		}
		// If there is no index in the node, we began the routine at an extremity
		// We found a border, let's keep going in the opposite way
		if( !( node->LineIndexes.size() ) )
		{
			*borderFound = !*borderFound;
		}
		else
		{
			end = true;
		}
	}

	// Else, we are in the middle of a line

	unsigned int i = 0;

	// If the node is not empty
	// Here, we decide which way to go through the line
	// By default, i = 0
	if( node->LineIndexes.size() )
	{	
		// If the current pixel is already in the container,
		// let i = 1 (opposite direction)
		if( node->LineIndexes.front() == it.GetIndex() || node->LineIndexes.back() == it.GetIndex() )
		{
			i = 1;
		}
		else
		{
			IndexType backIndex = node->LineIndexes.back();
			IndexType frontIndex = node->LineIndexes.front();

			// If the current pixel is not in the container, we check
			// if the default neighbor (i=0) is in the container.
			// If so, we go to the opposite direction (i=1)
			i = ( backIndex == it.GetIndex( map[0] ) || frontIndex == it.GetIndex( map[0] ) );

			// If opposite direction (i=1), and count > 1 (so in the middle of a line)
			// and if the opposite direction is already in the container, the
			// line is entirely scanned, we can quit.
			if( i && count > 1 && ( backIndex == it.GetIndex( map[1] ) || frontIndex == it.GetIndex( map[1] ) ) )
			{
				end = true;
				*borderFound = true;
			}
		}
	}

	// If the pixel where we are has a different identifier than the node
	// (we add +1 because the first node identifier is 0)
	// We can add the pixel to the container
	if( iit.GetCenterPixel() != node->Identifier + 1 )
	{
		// Depending of borderFound : push back or front
		if( *borderFound )
		{
			node->LineIndexes.push_back( it.GetIndex() );
			node->VesselnessValues.push_back( vit.GetCenterPixel() );
		}
		else
		{
			node->LineIndexes.push_front( it.GetIndex() );
			node->VesselnessValues.push_front( vit.GetCenterPixel() );
		}
		
		// We went here
		iit.SetCenterPixel( node->Identifier + 1 );
	}

	if( end )
	{
		return;
	}

	// We go to the next pixel
	it += it.GetOffset( map[i] );
	lit += lit.GetOffset( map[i] );
	iit += iit.GetOffset( map[i] );
	vit += vit.GetOffset( map[i] );

	this->MakeNode( it, lit, iit, vit, node, borderFound );
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
void SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::ConnectSelf( NodePointerType node )
{
	EdgePointerType edge = this->GetOutput()->CreateNewEdge( node, node, 1 );
	edge->ReverseEdgeIdentifier = edge->Identifier;
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
void SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::ConnectNodes( ImageIteratorType & it, LabelImageIteratorType & lit, IdentifierImageIteratorType & iit, VesselnessImageIteratorType & vit, NodePointerType node, LabelType label )
{
	ArrayType map;

	unsigned int count = this->GetNeighborhoodMap( it, map );

	// If we are not in an intersection
	if( count <= 2 )
	{
		// If the pixel is not virgin
		if( iit.GetCenterPixel() )
		{
			NodePointerType tmp = this->GetOutput()->GetNodePointer( iit.GetCenterPixel() - 1 );
			// We assign the label of the extremities
			// We need this label to know in which direction nodes
			// are connected with each other
			if( tmp->LineIndexes.size() < 2 || iit.GetIndex() == *( ++( tmp->LineIndexes.begin() ) ) )
			{
				tmp->FrontLabel = m_Label;
			}
			else
			{
				tmp->BackLabel = m_Label;
			}
		}
	}


	if( lit.GetCenterPixel() == label )
	{
		return;
	}

	// If we are not in an intersection
	if( count <= 2 )
	{
		// Not virgin
		if( iit.GetCenterPixel() )
		{
			bool nodeAlreadyProcessed = false;
			// We get the current node pointer
			NodePointerType nodetmp = this->GetOutput()->GetNodePointer( iit.GetCenterPixel()-1 );
			// For every incoming edges
			for( unsigned int i = 0; i < nodetmp->IncomingEdges.size(); ++i )
			{
				if( this->GetOutput()->GetSourceNodePointer( nodetmp->IncomingEdges[i] )->Identifier == node->Identifier
				|| this->GetOutput()->GetSourceNodePointer( nodetmp->IncomingEdges[i] )->Identifier == node->Identifier )
				{
					nodeAlreadyProcessed = true;
				}
			}
			for( unsigned int i = 0; i < nodetmp->OutgoingEdges.size(); ++i )
			{
				if( this->GetOutput()->GetSourceNodePointer( nodetmp->OutgoingEdges[i] )->Identifier == node->Identifier
				|| this->GetOutput()->GetSourceNodePointer( nodetmp->OutgoingEdges[i] )->Identifier == node->Identifier )
				{
					nodeAlreadyProcessed = true;
				}
			}
			
			if( !nodeAlreadyProcessed )
			{
				// We connect node and nodetmp
				EdgePointerType edge2 = this->GetOutput()->CreateNewEdge( nodetmp, node, 1 );
				EdgePointerType edge1 = this->GetOutput()->CreateNewEdge( node, nodetmp, 1 );
				edge1->ReverseEdgeIdentifier = edge2->Identifier;
				edge2->ReverseEdgeIdentifier = edge1->Identifier;
			}
		}
		return;
	}

	bool keepGoing = false;
	// We update the label of the intersection, that changed every time we pass by
	if( lit.GetCenterPixel() != label )
	{
		keepGoing = true;
		lit.SetCenterPixel( label );
	}


	for( unsigned int i = 0; i < count; ++i )
	{
		ImageIteratorType tit = it;
		LabelImageIteratorType tlit = lit;
		IdentifierImageIteratorType tiit = iit;
		VesselnessImageIteratorType tvit = vit;

		tit += tit.GetOffset(map[i]);
		tlit += tit.GetOffset(map[i]);
		tiit += tit.GetOffset(map[i]);
		tvit += tit.GetOffset(map[i]);


		if( keepGoing )
		{
			this->ConnectNodes( tit, tlit, tiit, tvit, node, label );
		}
	}

}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
void SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::FillIntersection( ImageIteratorType & it, LabelImageIteratorType & lit, IdentifierImageIteratorType & iit, VesselnessImageIteratorType & vit, LabelType & label )
{
	// If pixel is not virgin
	if( lit.GetCenterPixel() )
	{
		return;
	}

	// We passed here
	lit.SetCenterPixel( label );

	ArrayType map;

	unsigned int count = this->GetNeighborhoodMap( it, map );

	if( count <= 2 )
	{
		return;
	}

	for( unsigned int i = 0; i < count; ++i )
	{
		ImageIteratorType tit = it;
		LabelImageIteratorType tlit = lit;
		IdentifierImageIteratorType tiit = iit;
		VesselnessImageIteratorType tvit = vit;

		tit += tit.GetOffset(map[i]);
		tlit += tit.GetOffset(map[i]);
		tiit += tit.GetOffset(map[i]);
		tvit += tit.GetOffset(map[i]);

		// Keep going

		this->FillIntersection( tit, tlit, tiit, tvit, label );
	}
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
unsigned int SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::GetNeighborhoodMap( const ImageIteratorType & it, ArrayType & map ) const
{
	unsigned int count = 0;
	for( unsigned int i = 0; i < m_NeighborhoodCount; ++i )
	{
		if( m_Region.IsInside( it.GetIndex( i ) ) && it.GetPixel( i ) && it.GetIndex( i ) != it.GetIndex() )
		{
			map.push_back( i );
			++count;
		}
	}
	return count;
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
unsigned int SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::GetNeighborhoodCount( const ImageIteratorType & it ) const
{
	unsigned int count = 0;
	for( unsigned int i = 0; i < m_NeighborhoodCount; ++i )
	{
		if( m_Region.IsInside( it.GetIndex( i ) ) && it.GetPixel( i ) )
		{
			++count;
		}
	}
	return count == 2 || count == 3;
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
bool SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::ReachBranchSkeleton( ImageIteratorType & it, LabelImageIteratorType & lit, IdentifierImageIteratorType & iit, VesselnessImageIteratorType & vit ) const
{
	while( !it.IsAtEnd() && ( !this->GetNeighborhoodCount( it ) || ( iit.GetCenterPixel() || !it.GetCenterPixel() ) ) )
	{
		++it;
		++iit;
		++lit;
		++vit;
	}
	return !it.IsAtEnd();
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
void SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::SetVesselnessInput( const TVesselnessImage* vesselness )
{
	this->SetNthInput( 1, const_cast< TVesselnessImage* > ( vesselness ) ); 
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
const TVesselnessImage * SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::GetVesselnessInput() const
{
	return static_cast< const TVesselnessImage * > ( this->ProcessObject::GetInput( 1 ) );
}

} // namespace itk

#endif
