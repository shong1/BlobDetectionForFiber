#ifndef __itkSkeletonGraph_hxx
#define __itkSkeletonGraph_hxx

#include "itkSkeletonGraph.h"
#include "itkBresenhamLine.h"

#include <iterator>
#include <algorithm>

namespace itk
{

template< class TGraphTraits >
typename SkeletonGraph< TGraphTraits >::Pointer SkeletonGraph< TGraphTraits >
::Copy() const
{
	Pointer output = Self::New();
	for( typename NodeContainerType::ConstIterator it = this->GetNodeContainer()->Begin(); it != this->GetNodeContainer()->End(); ++it )
	{
		NodePointerType node = output->CreateNewNode();
		node->LineIndexes = it.Value().LineIndexes;
		node->VesselnessValues = it.Value().VesselnessValues;
		node->BackLabel = it.Value().BackLabel;
		node->FrontLabel = it.Value().FrontLabel;
		node->Weight = it.Value().Weight;
		node->Vector = it.Value().Vector;
		node->Point = it.Value().Point;
		node->MeanVesselness = it.Value().MeanVesselness;
		node->Force = it.Value().Force;
		node->AverageDistance = it.Value().AverageDistance;
		node->TheoricalAverageDistance = it.Value().TheoricalAverageDistance;
	}

	for( typename NodeContainerType::ConstIterator it = this->GetNodeContainer()->Begin(), outIt = output->GetNodeContainer()->Begin(); it != this->GetNodeContainer()->End(); ++it, ++outIt )
	{
		for( unsigned int i = 0; i < it.Value().OutgoingEdges.size(); ++i )
		{
			const NodeIdentifierType target = this->GetTargetNode( it.Value().OutgoingEdges[i] ).Identifier;
			EdgePointerType edge = output->CreateNewEdge( it.Value().Identifier, target );
			const EdgeType tmp = this->GetEdge( it.Value().Identifier, target );
			edge->Weight = tmp.Weight;
			edge->f = tmp.f;
			edge->Vector = tmp.Vector;
			edge->Point = tmp.Point;
			edge->AverageDistance = tmp.AverageDistance;
			edge->TheoricalAverageDistance = tmp.TheoricalAverageDistance;
			edge->LineIndexes = tmp.LineIndexes;
		}
	}

	for( typename NodeContainerType::Iterator it = output->GetNodeContainer()->Begin(); it != output->GetNodeContainer()->End(); ++it )
	{
		output->SynchronizeIncomingOnOutgoingEdges( & ( it.Value() ) );
	}

	output->ConnectNodeExtremities();

	return output;
}

template< class TGraphTraits >
void SkeletonGraph< TGraphTraits >
::InitializeGraphAttributes()
{
	for( typename NodeContainerType::Iterator it = this->GetNodeContainer()->Begin(); it != this->GetNodeContainer()->End(); ++it )
	{
		this->InitializeNode( & ( it.Value() ) );
	}
	for( typename NodeContainerType::Iterator it = this->GetNodeContainer()->Begin(); it != this->GetNodeContainer()->End(); ++it )
	{
		for( unsigned int j = 1; j < it.Value().OutgoingEdges.size(); ++j )
		{
			this->InitializeEdge( this->GetEdgePointer( it.Value().OutgoingEdges[j] ), it.Value(), this->GetTargetNode( it.Value().OutgoingEdges[j] ) );
		}
    }    
	for( typename NodeContainerType::Iterator it = this->GetNodeContainer()->Begin(); it != this->GetNodeContainer()->End(); ++it )
	{
		this->InitializeNodeInteraction( & ( it.Value() ) );
		this->InitializeOutgoingEdgesInteraction( it.Value() );
	}
}

template< class TGraphTraits >
void SkeletonGraph< TGraphTraits >
::InitializeOutgoingEdgesInteraction( const NodeType & node )
{
	this->SetOutgoingEdgesWeight( node );
}

template< class TGraphTraits >
void SkeletonGraph< TGraphTraits >
::InitializeNodeInteraction( NodePointerType node )
{
	node->Force = this->Force( *node );
	node->Weight = this->Weight( *node );
}

template< class TGraphTraits >
void SkeletonGraph< TGraphTraits >
::InitializeEdge( EdgePointerType edge, const NodeType & source, const NodeType & target )
{
	edge->LineIndexes = this->MergeIndexes( source, target, this->IsConnectedToFront( source, target ), this->IsConnectedToFront( target, source ) );

	edge->Point = this->P< typename IndexContainerType::const_iterator >( edge->LineIndexes.begin(), edge->LineIndexes.end() );
	edge->Vector = this->V< typename IndexContainerType::const_iterator >( edge->LineIndexes.begin(), edge->LineIndexes.end(), edge->Point );
	edge->AverageDistance = this->AverageDistance< typename IndexContainerType::const_iterator >( edge->LineIndexes.begin(), edge->LineIndexes.end(), edge->Point, edge->Vector );
	edge->TheoricalAverageDistance = this->TheoricalAverageDistance( edge->Vector, 2 * edge->LineIndexes.size() );
	edge->f = this->f( source, target );
}

template< class TGraphTraits >
void SkeletonGraph< TGraphTraits >
::InitializeNode( NodePointerType node )
{
	node->Point = this->P< typename IndexContainerType::const_iterator >( node->LineIndexes.begin(), node->LineIndexes.end() );
	node->Vector = this->V< typename IndexContainerType::const_iterator >( node->LineIndexes.begin(), node->LineIndexes.end(), node->Point );
	node->AverageDistance = this->AverageDistance< typename IndexContainerType::const_iterator >( node->LineIndexes.begin(), node->LineIndexes.end(), node->Point, node->Vector );
	node->TheoricalAverageDistance = this->TheoricalAverageDistance( node->Vector, 2 * node->LineIndexes.size() );
	node->MeanVesselness = this->Mean< typename ScalarContainerType::const_iterator >( node->VesselnessValues.begin(), node->VesselnessValues.end() );

}

template< class TGraphTraits >
void SkeletonGraph< TGraphTraits >
::ConnectNodeExtremities()
{
	for( typename NodeContainerType::Iterator it = this->GetNodeContainer()->Begin(); it != this->GetNodeContainer()->End(); ++it )
	{
		if( it.Value().OutgoingEdges.size() > 1 )
		{
			for( unsigned int i = 1; i < it.Value().OutgoingEdges.size(); ++i )
			{
				NodeType target = this->GetTargetNode( it.Value().OutgoingEdges[i] );
				if( this->IsConnectedToFront( it.Value(), target ) )
				{
					it.Value().OutgoingEdgesFront.push_back( it.Value().OutgoingEdges[i] );
					it.Value().IncomingEdgesFront.push_back( it.Value().IncomingEdges[i] );
				}
				if( this->IsConnectedToBack( it.Value(), target ) )
				{
					it.Value().OutgoingEdgesBack.push_back( it.Value().OutgoingEdges[i] );
					it.Value().IncomingEdgesBack.push_back( it.Value().IncomingEdges[i] );
				}
			}
			it.Value().OutgoingEdges.erase( it.Value().OutgoingEdges.begin() + 1, it.Value().OutgoingEdges.end() );
			it.Value().IncomingEdges.erase( it.Value().IncomingEdges.begin() + 1, it.Value().IncomingEdges.end() );
			std::copy( it.Value().OutgoingEdgesFront.begin(), it.Value().OutgoingEdgesFront.end(), std::back_inserter( it.Value().OutgoingEdges ) );
			std::copy( it.Value().IncomingEdgesFront.begin(), it.Value().IncomingEdgesFront.end(), std::back_inserter( it.Value().IncomingEdges ) );
			std::copy( it.Value().OutgoingEdgesBack.begin(), it.Value().OutgoingEdgesBack.end(), std::back_inserter( it.Value().OutgoingEdges ) );
			std::copy( it.Value().IncomingEdgesBack.begin(), it.Value().IncomingEdgesBack.end(), std::back_inserter( it.Value().IncomingEdges ) );
		}
	}
}

template< class TGraphTraits >
void SkeletonGraph< TGraphTraits >
::SetOutgoingEdgesWeight( const NodeType & node )
{
	double denom = 0;
	for( unsigned int j = 0; j < node.OutgoingEdges.size(); ++j )
	{
		NodeType target = this->GetTargetNode( node.OutgoingEdges[j] );
		denom += target.Force / node.MeanVesselness * target.MeanVesselness;
	}
	for( unsigned int j = 0; j < node.OutgoingEdges.size(); ++j )
	{
		NodeType target = this->GetTargetNode( node.OutgoingEdges[j] );
		this->GetEdge( node.OutgoingEdges[j] ).Weight = target.Force * node.MeanVesselness / target.MeanVesselness / denom;
	}
}

template< class TGraphTraits >
typename SkeletonGraph< TGraphTraits >
::NodePointerType SkeletonGraph< TGraphTraits >
::GetMaximumNodeWeight( const NodeWeightType & threshold )
{
	double max = 0;
	NodePointerType node;
	bool maximumFound = false;
	for( typename NodeContainerType::Iterator it = this->GetNodeContainer()->Begin(); it != this->GetNodeContainer()->End(); ++it )
	{
		if( it.Value().Weight > max && ( threshold < 0 || threshold > it.Value().Weight ) )
		{
			node = & ( it.Value() );
			maximumFound = true;
			max = it.Value().Weight;
		}
	}
	if( !maximumFound )
	{
		return NULL;
	}
	return node;
}

template< class TGraphTraits >
typename SkeletonGraph< TGraphTraits >
::NodeWeightType SkeletonGraph< TGraphTraits >
::Weight( const NodeType & node ) const
{
	return node.Force * node.MeanVesselness;
}

template< class TGraphTraits >
double SkeletonGraph< TGraphTraits >
::f( const NodeType & i, const NodeType & j ) const
{
	// I and J serve to know if the two nodes have an
	// opposite direction or not
    double cos = std::abs( i.Vector * j.Vector );
    VectorType I, J;
    if( this->IsConnectedToFront( i, j ) )
    {
		for( unsigned int k = 0; k < Dimension; ++k )
		{
			I[k] = i.LineIndexes.front()[k];
		}
    }
    else
    {
		for( unsigned int k = 0; k < Dimension; ++k )
		{
			I[k] = i.LineIndexes.back()[k];
		}
    }
    if( this->IsConnectedToFront( j, i ) )
    {
		for( unsigned int k = 0; k < Dimension; ++k )
		{
			J[k] = j.LineIndexes.front()[k];
		}
    }
    else
    {
		for( unsigned int k = 0; k < Dimension; ++k )
		{
			J[k] = j.LineIndexes.back()[k];
		}
    }
	double theta = std::acos( cos );
	if( theta != theta )
	{
		theta = 0;
	}

	double result = std::exp( -2 * theta * std::abs( std::tan( theta / 2 ) ) );
	if( result != result )
	{
		return 0;
	}
    if( ( i.Point - I ) * ( j.Point - J ) > 0 )
    {
        result = -result;
    }
	return result;
}

template< class TGraphTraits >
typename SkeletonGraph< TGraphTraits >
::IndexContainerType SkeletonGraph< TGraphTraits >
::MergeIndexes( const NodeType & i, const NodeType & j, bool ifront, bool jfront ) const
{
	IndexContainerType indexes;
	std::copy( i.LineIndexes.begin(), i.LineIndexes.end(), std::back_inserter( indexes ) );
	if( ifront && jfront )
	{
		if( i.LineIndexes.front() == j.LineIndexes.front() )
		{
			std::copy( ++( j.LineIndexes.begin() ), j.LineIndexes.end(), std::front_inserter( indexes ) );
		}
		else
		{
			std::copy( j.LineIndexes.begin(), j.LineIndexes.end(), std::front_inserter( indexes ) );
		}
	}
	if( ifront && !jfront )
	{
		if( i.LineIndexes.front() == j.LineIndexes.back() )
		{
				std::copy( ++( j.LineIndexes.rbegin() ), j.LineIndexes.rend(), std::front_inserter( indexes ) );
		}
		else
		{
			std::copy( j.LineIndexes.rbegin(), j.LineIndexes.rend(), std::front_inserter( indexes ) );
		}
	}
	if( !ifront && jfront )
	{
		if( i.LineIndexes.back() == j.LineIndexes.front() )
		{
			std::copy( ++( j.LineIndexes.begin() ), j.LineIndexes.end(), std::back_inserter( indexes ) );
		}
		else
		{
			std::copy( j.LineIndexes.begin(), j.LineIndexes.end(), std::back_inserter( indexes ) );
		}
	}
	if( !ifront && !jfront )
	{
		if( i.LineIndexes.back() == j.LineIndexes.back() )
		{
			std::copy( ++( j.LineIndexes.rbegin() ), j.LineIndexes.rend(), std::back_inserter( indexes ) );
		}
		else
		{
			std::copy( j.LineIndexes.rbegin(), j.LineIndexes.rend(), std::back_inserter( indexes ) );
		}
	}
	return indexes;
}

template< class TGraphTraits >
template< class TIterator >
typename SkeletonGraph< TGraphTraits >
::VectorType SkeletonGraph< TGraphTraits >
::P( TIterator begin, TIterator end ) const
{
	VectorType p;
	p.Fill( 0 );
	for( TIterator it = begin; it != end; ++it )
	{
		for( unsigned int i = 0; i < Dimension; ++i )
		{
			p[i] += (*it)[i];
		}
	}
	return  p / ( double ) std::distance( begin, end );
}

// if M = B + t * u (u is a vector) defines a line
// the distance of the point A to the line is :
// d = ||A-B-t*u||
// t = (BA.u) / ||u||^2 
// . is the inner product
template< class TGraphTraits >
template< class TIterator >
typename SkeletonGraph< TGraphTraits >
::VectorType SkeletonGraph< TGraphTraits >
::V( TIterator begin, TIterator end, const VectorType & p ) const
{
	VectorType v;
	v.Fill( 0 );
	VectorType d;
	d.Fill( 0 );

	VectorType tmp;

	unsigned int count = 0;

	for( TIterator it = begin; it != end; ++it )
	{
		for( unsigned int i = 0; i < Dimension; ++i )
		{
			tmp[i] = (*it)[i] - p[i];
			v[i] += tmp[i] * tmp[i];
		}
		++count;
		if( count > std::distance( begin, end ) / 2 )
		{
			for( unsigned int i = 0; i < Dimension; ++i )
			{
				d[i] += (*it)[i];
			}
		}
	}

	d /= ( count + 1 ) / 2;
	d -= p;

	for( unsigned int i = 0; i < Dimension; ++i )
	{
		v[i] /= std::distance( begin, end );
		v[i] = std::sqrt( v[i] );
		if( d[i] < 0 )
		{
			v[i] = -v[i];
		}
	}
	v.Normalize();
	return v;
}

template< class TGraphTraits >
template< class TIterator >
double SkeletonGraph< TGraphTraits >
::AverageDistance( TIterator begin, TIterator end, const VectorType & p, const VectorType & v ) const
{
	double dist = 0;
	for( TIterator it = begin; it != end; ++it )
	{
		dist += this->Distance( *it, p, v );
	}
	return dist / ( double ) std::distance( begin, end );
}

template< class TGraphTraits >
double SkeletonGraph< TGraphTraits >
::Distance( const IndexType & index, const VectorType & p, const VectorType & v ) const
{
	VectorType A;
	float t = 0;
	for( unsigned int i = 0; i < Dimension; ++i ) 
	{
		A[i] = index[i];
		t += v[i] * ( p[i] - A[i] );
	}

	VectorType M = p - t * v;

	return ( A - M ).GetNorm();
}

template< class TGraphTraits >
double SkeletonGraph< TGraphTraits >
::TheoricalAverageDistance( const VectorType & v, double norm ) const
{
	typedef BresenhamLine< Dimension > LineType;
	LineType line;
	typename LineType::IndexType p0, p1;
	VectorType ptmp;
	p0.Fill( 0 );
	for( unsigned int i = 0; i < Dimension; ++i )
	{
		p1[i] = v[i] * norm;
	}
	typename LineType::IndexArray indexes = line.BuildLine( p0, p1 );
	
	return this->AverageDistance< typename LineType::IndexArray::const_iterator >( indexes.begin(), indexes.end(), this->P< typename LineType::IndexArray::const_iterator >( indexes.begin(), indexes.end() ), v );
}

template< class TGraphTraits >
template< class TIterator >
double SkeletonGraph< TGraphTraits >
::Mean( TIterator begin, TIterator end ) const
{
	double value = 0;
	for( TIterator it = begin; it != end; ++it )
	{
		value += *it;
	}
	return value / ( double ) std::distance( begin, end );
}

template< class TGraphTraits >
double SkeletonGraph< TGraphTraits >
::Force( const NodeType & node ) const
{
	double max = 0;
	unsigned int j = 1;
	for( j = 1; j < node.OutgoingEdges.size(); ++j )
	{
		double tmp = this->Force( node, this->GetTargetNode( node.OutgoingEdges[j] ), this->GetEdge( node.OutgoingEdges[j] ) );
		max = tmp > max ? tmp : max;
	}
	if( j == 1 )
	{
		return 0;
	}
	return max;
}

template< class TGraphTraits >
double SkeletonGraph< TGraphTraits >
::Force( const NodeType & i, const NodeType & j, const EdgeType & ij ) const
{
	double x = 1;
	for( unsigned int k = 0; k < 2; ++k )
	{
		x *= ( 1 - ij.f ) + std::abs( 1 - std::log( 1 + ( 1 + std::abs( ij.AverageDistance - ij.TheoricalAverageDistance ) )  / ( 1 + std::abs( i.AverageDistance - i.TheoricalAverageDistance ) ) ) );
	}

	double y = ( double ) i.LineIndexes.size() / ( double ) j.LineIndexes.size();
	double D = std::abs( ij.AverageDistance - ij.TheoricalAverageDistance );
	D = D * D;
	y = 0;
	if( y < 1 )
	{
		y = 0;
	}
	else
	{
		y = std::exp( -1 / ( y - 1 ) * ( y - 1 ) );
	}
	y = y > ij.f ? y : ij.f;

	return y / ( 1 + D );
}

template< class TGraphTraits >
bool SkeletonGraph< TGraphTraits >
::IsConnectedToFront( const NodeType & i, const NodeType & j ) const
{
	return i.FrontLabel && ( i.FrontLabel == j.FrontLabel || i.FrontLabel == j.BackLabel );
}

template< class TGraphTraits >
bool SkeletonGraph< TGraphTraits >
::IsConnectedToBack( const NodeType & i, const NodeType & j ) const
{
	return i.BackLabel && ( i.BackLabel == j.FrontLabel || i.BackLabel == j.BackLabel );
}

template< class TGraphTraits >
bool SkeletonGraph< TGraphTraits >
::IsOnlyConnectedToFront( const NodeType & i, const NodeType & j ) const
{
	return this->IsConnectedToFront( i, j ) && ( !this->IsConnectedToBack( i, j ) || !i.BackLabel );
}

template< class TGraphTraits >
bool SkeletonGraph< TGraphTraits >
::IsOnlyConnectedToBack( const NodeType & i, const NodeType & j ) const
{
	return ( !this->IsConnectedToFront( i, j ) || !i.FrontLabel ) && this->IsConnectedToBack( i, j );
}

template< class TGraphTraits >
bool SkeletonGraph< TGraphTraits >
::IsSelfConnected( const NodeType & i ) const
{
	return i.BackLabel == i.FrontLabel && i.BackLabel || !i.BackLabel && !i.FrontLabel;
}

template< class TGraphTraits >
void SkeletonGraph< TGraphTraits >
::SynchronizeIncomingOnOutgoingEdges( NodePointerType node )
{
	EdgeIdentifierContainerType buf;
	for( unsigned int j = 0; j < node->OutgoingEdges.size(); ++j )
	{
		NodeType tmp = this->GetTargetNode( node->OutgoingEdges[j] );
		for( unsigned int i = 0; i < node->IncomingEdges.size(); ++i )
		{
			if( this->GetSourceNode( node->IncomingEdges[i] ).Identifier == tmp.Identifier )
			{
				buf.push_back( node->IncomingEdges[i] );
			}
		}
	}
	node->IncomingEdges.swap( buf );
	buf.clear();
}

template< class TGraphTraits >
void SkeletonGraph< TGraphTraits >
::DeleteLinks( NodePointerType & i, bool front )
{
	EdgeIdentifierContainerType *ijOut, *ijIn;
	if( front )
	{
		ijOut = &(i->OutgoingEdgesFront);
		ijIn = &(i->IncomingEdgesFront);
	}
	else
	{
		ijOut = &(i->OutgoingEdgesBack);
		ijIn = &(i->IncomingEdgesBack);
	}
	for( unsigned int k = 0; k < ijOut->size(); ++k )
	{
		NodePointerType j = this->GetTargetNodePointer( (*ijOut)[k] );
		bool nodeOutgoingSide = this->IsConnectedToFront( *j, *i );
		if( nodeOutgoingSide && !this->IsOnlyConnectedToFront( *j, *i ) )
		{
			if( j->FrontLabel && ( front && i->FrontLabel == j->FrontLabel || !front && i->BackLabel == j->FrontLabel ) )
			{
				nodeOutgoingSide = true;
			}
			if( j->BackLabel && ( front && i->FrontLabel == j->BackLabel || !front && i->BackLabel == j->BackLabel ) )
			{
				nodeOutgoingSide = false;
			}
		}
		unsigned int begin, end;
		if( nodeOutgoingSide )
		{
			begin = 1;
			end = 1 + j->OutgoingEdgesFront.size();
		}
		else
		{
			begin = 1 + j->OutgoingEdgesFront.size();
			end = j->OutgoingEdges.size();
		}
		for( unsigned int l = begin; l < end; ++l )
		{
			if( this->GetTargetNode( j->OutgoingEdges[l] ).Identifier == i->Identifier )
			{
				j->OutgoingEdges.erase( j->OutgoingEdges.begin() + l );
				j->IncomingEdges.erase( j->IncomingEdges.begin() + l );
				--l;
				--end;
			}
		}
		EdgeIdentifierContainerType *outgoing, *incoming;
		if( nodeOutgoingSide )
		{
			outgoing = &(j->OutgoingEdgesFront);
			incoming = &(j->IncomingEdgesFront);
		}
		else
		{
			outgoing = &(j->OutgoingEdgesBack);
			incoming = &(j->IncomingEdgesBack);
		}
		for( unsigned int l = 0; l < outgoing->size(); ++l )
		{
			if( this->GetTargetNode( (*outgoing)[l] ).Identifier == i->Identifier )
			{
				outgoing->erase( outgoing->begin() + l );
				incoming->erase( incoming->begin() + l );
			}
		}
	}
	if( front )
	{
		i->OutgoingEdges.erase( i->OutgoingEdges.begin() + 1, i->OutgoingEdges.begin() + 1 + ijOut->size() );
		i->IncomingEdges.erase( i->IncomingEdges.begin() + 1, i->IncomingEdges.begin() + 1 + ijOut->size() );
	}
	else
	{
		i->OutgoingEdges.erase( i->OutgoingEdges.begin() + i->OutgoingEdgesFront.size() + 1, i->OutgoingEdges.end() );
		i->IncomingEdges.erase( i->IncomingEdges.begin() + i->IncomingEdgesFront.size() + 1, i->IncomingEdges.end() );
	}
	ijOut->clear();
	ijIn->clear();
}

template< class TGraphTraits >
double SkeletonGraph< TGraphTraits >
::ConnectionLikelihood( const NodeType & i, const NodeType & j, const EdgeType & edge ) const
{
	return edge.f;
}

template< class TGraphTraits >
double SkeletonGraph< TGraphTraits >
::MaxDistance( const IndexContainerType & indexes, const VectorType & point, const VectorType & vector, unsigned int & pos ) const
{
	double dist = 0;
	unsigned int i = 0;
	for( typename IndexContainerType::const_iterator it = indexes.begin(); it != indexes.end(); ++it, ++i )
	{
		double tmp = this->Distance( *it, point, vector );
		if( tmp > dist )
		{
			dist = tmp;
			pos = i;
		}
	}
	return dist;
}

template< class TGraphTraits >
bool SkeletonGraph< TGraphTraits >
::MergeSuccess( NodePointerType node, const double & threshold )
{
	if( node == NULL || this->IsSelfConnected( *node ) )
	{
		return false;
	}

	double connectionLikelihoodMax = 0;
	unsigned int iMax = 0;
	bool iFront, jFront = true;
	bool iFound = false;


	for( unsigned int i = 0; i < node->OutgoingEdgesFront.size(); ++i )
	{
		double connectionLikelihood = this->ConnectionLikelihood( *node, this->GetTargetNode( node->OutgoingEdgesFront[i] ), this->GetEdge( node->OutgoingEdgesFront[i] ) );
		if( connectionLikelihoodMax < connectionLikelihood )
		{
			connectionLikelihoodMax = connectionLikelihood;
			iMax = i;
			iFront = true;
			iFound = true;
		}
	}
	for( unsigned int i = 0; i < node->OutgoingEdgesBack.size(); ++i )
	{
		double connectionLikelihood = this->ConnectionLikelihood( *node, this->GetTargetNode( node->OutgoingEdgesBack[i] ), this->GetEdge( node->OutgoingEdgesBack[i] ) );
		if( connectionLikelihoodMax < connectionLikelihood )
		{
			connectionLikelihoodMax = connectionLikelihood;
			iMax = i;
			iFront = false;
			iFound = true;
		}
	}
	if( connectionLikelihoodMax < 0.5 )
	{
		iFound = false;
	}
	if( iFound )
	{
		NodePointerType target;
		if( iFront )
		{
			target = this->GetTargetNodePointer( node->OutgoingEdgesFront[iMax] );
		}
		else
		{
			target = this->GetTargetNodePointer( node->OutgoingEdgesBack[iMax] );
		}
		jFront = this->IsConnectedToFront( *target, *node );

		if( this->CanBeMerged( *node, *target, threshold ) )
		{
			if( this->Merge( node, target, iFront, jFront ) )
			{
				this->UpdateNode( node );
				return true;
			}
		}
	}
	return false;
}

template< class TGraph >
void SkeletonGraph< TGraph >
::UpdateNode( NodePointerType node )
{
	this->InitializeNode( node );
	for( unsigned int j = 1; j < node->OutgoingEdges.size(); ++j )
	{
		this->InitializeEdge( this->GetEdgePointer( node->OutgoingEdges[j] ), *node, this->GetTargetNode( node->OutgoingEdges[j] ) );
	}
	for( unsigned int j = 1; j < node->OutgoingEdges.size(); ++j )
	{
		NodePointerType target = this->GetTargetNodePointer( node->OutgoingEdges[j] );
		this->InitializeNodeInteraction( target );
		this->InitializeOutgoingEdgesInteraction( *target );
	}
}


template< class TGraph >
bool SkeletonGraph< TGraph >
::CanBeMerged( const NodeType & source, const NodeType & target, double threshold ) const
{
	IndexContainerType mergedIndexes = this->MergeIndexes( source, target, this->IsConnectedToFront( source, target ), this->IsConnectedToFront( target, source ) );
	VectorType point = this->P< typename IndexContainerType::const_iterator >( mergedIndexes.begin(), mergedIndexes.end() );
	VectorType vector = this->V< typename IndexContainerType::const_iterator >( mergedIndexes.begin(), mergedIndexes.end(), point );
	unsigned int junk = 0;
	double max = this->MaxDistance( mergedIndexes, point, vector, junk );
	mergedIndexes.clear();
	return max < threshold;
}

template< class TGraphTraits >
bool SkeletonGraph< TGraphTraits >
::IsLikelyConnectedToFront( const NodeType & i, const NodeType & j ) const
{
	IndexType ifront = i.LineIndexes.front(), iback = i.LineIndexes.back(), jdx = j.LineIndexes.front();
	unsigned int frontdist = 0, backdist = 0;
	for( unsigned int k = 0; k < Dimension; ++k )
	{
		frontdist += ( ifront[k] - jdx[k] ) * ( ifront[k] - jdx[k] );
		backdist += ( iback[k] - jdx[k] ) * ( iback[k] - jdx[k] );
	}
	return frontdist < backdist;
}

template< class TGraphTraits >
bool SkeletonGraph< TGraphTraits >
::IsLikelyConnectedToBack( const NodeType & i, const NodeType & j ) const
{
	return !this->IsLikelyConnetedToFront( i, j );
}

template< class TGraphTraits >
double SkeletonGraph< TGraphTraits >
::VirtualConnectionLikelihood( const NodeType & i, const NodeType & j ) const
{
	return 1;
}

template< class TGraphTraits >
void SkeletonGraph< TGraphTraits >
::FillGaps( double threshold )
{
	for( typename NodeContainerType::Iterator it = this->GetNodeContainer()->Begin(); it != this->GetNodeContainer()->End(); ++it )
	{
		double likelihood = 0;
		NodePointerType j;
		if( it.Value().LineIndexes.size() )
		{
			typename NodeContainerType::Iterator jt = it;
			unsigned int count = 0;
			for( ++jt ; jt != this->GetNodeContainer()->End(); ++jt )
			{
				++count;
				if( jt.Value().Identifier != it.Value().Identifier && jt.Value().LineIndexes.size() && !this->AreConnected( it.Value(), jt.Value() ) )
				{
					double tmp = this->VirtualConnectionLikelihood( it.Value(), jt.Value() );
					if( tmp > likelihood )
					{
						likelihood = tmp;
						j = &( jt.Value() );
					}
				}
			}
			if( likelihood > threshold )
			{
				NodePointerType i = &( it.Value() );
				this->Connect( i, j, this->IsLikelyConnectedToFront( it.Value(), *j ), this->IsLikelyConnectedToFront( *j, it.Value() ) );
			}
		}
	}
}

template< class TGraphTraits >
bool SkeletonGraph< TGraphTraits >
::AreConnected( const NodeType & i, const NodeType &  j ) const
{
	for( unsigned int k = 1; k < i.OutgoingEdges.size(); ++k )
	{
		if( this->GetTargetNode( i.OutgoingEdges[k] ).Identifier == j.Identifier )
		{
			return true;
		}	
	}
	return false;
}

template< class TGraphTraits >
void SkeletonGraph< TGraphTraits >
::Connect( NodePointerType & i, NodePointerType & j, bool iFront, bool jFront )
{
	unsigned int ibegin, iend, jbegin, jend;
	EdgePointerType ijEdge = this->CreateNewEdge( i->Identifier, j->Identifier );
	EdgePointerType jiEdge = this->CreateNewEdge( j->Identifier, i->Identifier );
	this->InitializeEdge( jiEdge, *j, *i );
	if( iFront )
	{
		i->OutgoingEdgesFront.push_back( ijEdge->Identifier );
		i->IncomingEdgesFront.push_back( jiEdge->Identifier );
		std::move( --( i->OutgoingEdges.end() ), --( i->OutgoingEdges.end() ), i->OutgoingEdges.begin() + 1 + i->OutgoingEdgesFront.size() );
		std::move( --( i->IncomingEdges.end() ), --( i->IncomingEdges.end() ), i->IncomingEdges.begin() + 1 + i->IncomingEdgesFront.size() );
	}
	else
	{
		i->OutgoingEdgesBack.push_back( ijEdge->Identifier );
		i->IncomingEdgesBack.push_back( jiEdge->Identifier );
	}
	if( jFront )
	{
		j->OutgoingEdgesFront.push_back( jiEdge->Identifier );
		j->IncomingEdgesFront.push_back( ijEdge->Identifier );
		std::move( --( j->OutgoingEdges.end() ), --( j->OutgoingEdges.end() ), j->OutgoingEdges.begin() + 1 + j->OutgoingEdgesFront.size() );
		std::move( --( j->IncomingEdges.end() ), --( j->IncomingEdges.end() ), j->IncomingEdges.begin() + 1 + j->IncomingEdgesFront.size() );
	}
	else
	{
		j->OutgoingEdgesBack.push_back( jiEdge->Identifier );
		j->IncomingEdgesBack.push_back( ijEdge->Identifier );
	}
	this->InitializeNode( i );
	this->InitializeNode( j );
	this->InitializeEdge( ijEdge, *i, *j );
	this->InitializeNodeInteraction( i );
	this->InitializeNodeInteraction( j );
	if( !this->Merge( i, j, iFront, jFront, false ) )
	{
		std::cout << " Big problem " << std::endl;
	}
}


template< class TGraphTraits >
bool SkeletonGraph< TGraphTraits >
::Merge( NodePointerType & i, NodePointerType & j, bool ijOutgoingSide, bool jiOutgoingSide, bool reloopingPossible )
{
	if( this->IsSelfConnected( *i ) )
	{
		i->BackLabel = 0;
		i->FrontLabel = 0;
		this->DeleteLinks( i, true );
		i->OutgoingEdgesFront.clear();
		i->IncomingEdgesFront.clear();
		i->OutgoingEdgesBack.clear();
		i->IncomingEdgesBack.clear();
		i->OutgoingEdges.erase( i->OutgoingEdges.begin() + 1, i->OutgoingEdges.end() );
		i->IncomingEdges.erase( i->IncomingEdges.begin() + 1, i->IncomingEdges.end() );
		return false;
	}
	if( this->IsSelfConnected( *j ) )
	{
		j->BackLabel = 0;
		j->FrontLabel = 0;
		this->DeleteLinks( j, true );
		j->OutgoingEdgesFront.clear();
		j->IncomingEdgesFront.clear();
		j->OutgoingEdgesBack.clear();
		j->IncomingEdgesBack.clear();
		j->OutgoingEdges.erase( j->OutgoingEdges.begin() + 1, j->OutgoingEdges.end() );
		j->IncomingEdges.erase( j->IncomingEdges.begin() + 1, j->IncomingEdges.end() );
		return false;
	}

	bool keepGoing = true;

	if( reloopingPossible && ijOutgoingSide && this->IsConnectedToBack( *i, *j ) )
	{
		if( i->FrontLabel == j->FrontLabel )
		{
			jiOutgoingSide = false;
		}
		this->DeleteLinks( i, false );
		this->DeleteLinks( j, !jiOutgoingSide );
		keepGoing = false;
	}

	EdgeIdentifierContainerType *ij = &(i->OutgoingEdgesFront), *ji = &(j->OutgoingEdgesFront);
	EdgeIdentifierContainerType *ijReverse = &(i->OutgoingEdgesBack), *jiReverse = &(j->OutgoingEdgesBack);

	if( !ijOutgoingSide )
	{
		ij = &(i->OutgoingEdgesBack);
		ijReverse = &(i->OutgoingEdgesFront);
	}
	if( !jiOutgoingSide )
	{
		ji = &(j->OutgoingEdgesBack);
		jiReverse = &(j->OutgoingEdgesFront);
	}

	EdgeIdentifierContainerType *ijIncoming = &(i->IncomingEdgesFront), *jiIncoming = &(j->IncomingEdgesFront);
	EdgeIdentifierContainerType *ijIncomingReverse = &(i->IncomingEdgesBack), *jiIncomingReverse = &(j->IncomingEdgesBack);
	if( !ijOutgoingSide )
	{
		ijIncoming = &(i->IncomingEdgesBack);
		ijIncomingReverse = &(i->IncomingEdgesFront);
	}
	if( !jiOutgoingSide )
	{
		jiIncoming = &(j->IncomingEdgesBack);
		jiIncomingReverse = &(j->IncomingEdgesFront);
	}

	this->DeleteLinks( i, ijOutgoingSide );
	this->DeleteLinks( j, jiOutgoingSide );

	ij->clear();
	ijIncoming->clear();

	for( unsigned int k = 0; k < jiReverse->size(); ++k )
	{
		this->GetEdge( (*jiReverse)[k] ).SourceIdentifier = i->Identifier;
		this->GetEdge( (*jiIncomingReverse)[k] ).TargetIdentifier = i->Identifier;
	}

	std::copy( jiReverse->begin(), jiReverse->end(), std::back_inserter( *ij ) );
	std::copy( jiIncomingReverse->begin(), jiIncomingReverse->end(), std::back_inserter( *ijIncoming ) );

	if( ijOutgoingSide )
	{
		std::copy( jiReverse->begin(), jiReverse->end(), std::inserter( i->OutgoingEdges, ++( i->OutgoingEdges.begin() ) ) );
		std::copy( jiIncomingReverse->begin(), jiIncomingReverse->end(), std::inserter( i->IncomingEdges, ++( i->IncomingEdges.begin() ) ) );
	}
	else
	{
		std::copy( jiReverse->begin(), jiReverse->end(), std::back_inserter( i->OutgoingEdges ) );
		std::copy( jiIncomingReverse->begin(), jiIncomingReverse->end(), std::back_inserter( i->IncomingEdges ) );
	}

	j->OutgoingEdges.clear();
	j->OutgoingEdgesBack.clear();
	j->OutgoingEdgesFront.clear();

	if( ijOutgoingSide && jiOutgoingSide )
	{
		i->FrontLabel = j->BackLabel;
	}
	if( !ijOutgoingSide && jiOutgoingSide )
	{
		i->BackLabel = j->BackLabel;
	}
	if( ijOutgoingSide && !jiOutgoingSide )
	{
		i->FrontLabel = j->FrontLabel;
	}
	if( !ijOutgoingSide && !jiOutgoingSide )
	{
		i->BackLabel = j->FrontLabel;
	}

	if( !jiOutgoingSide && !ijOutgoingSide || jiOutgoingSide && ijOutgoingSide )
	{
		j->LineIndexes.reverse();
		j->VesselnessValues.reverse();
	}
	if( ijOutgoingSide )
	{
		i->LineIndexes.splice( i->LineIndexes.begin(), j->LineIndexes, j->LineIndexes.begin(), j->LineIndexes.end() );
		i->VesselnessValues.splice( i->VesselnessValues.end(), j->VesselnessValues, j->VesselnessValues.begin(), j->VesselnessValues.end() );
	}
	else
	{
		i->LineIndexes.splice( i->LineIndexes.end(), j->LineIndexes, j->LineIndexes.begin(), j->LineIndexes.end() );
		i->VesselnessValues.splice( i->VesselnessValues.end(), j->VesselnessValues, j->VesselnessValues.begin(), j->VesselnessValues.end() );
	}

	i->LineIndexes.unique();

	return keepGoing;
}

template< class TGraphTraits >
typename SkeletonGraph< TGraphTraits >
::NodeIdentifierType SkeletonGraph< TGraphTraits >
::MaxIntersectionLabel() const
{
	NodeIdentifierType id = 0;
	for( typename NodeContainerType::ConstIterator it = this->GetNodeContainer()->Begin(); it != this->GetNodeContainer()->End(); ++it )
	{
		if( id < it.Value().FrontLabel )
		{
			id = it.Value().FrontLabel;
		}
		if( id < it.Value().BackLabel )
		{
			id = it.Value().BackLabel;
		}
	}
	return id;
}

template< class TGraphTraits >
bool SkeletonGraph< TGraphTraits >
::CanBeSplit( const IndexContainerType & indexes, const VectorType & point, const VectorType & vector, double threshold, unsigned int & pos ) const
{
	double dist = threshold;
	unsigned int i = 0;
	unsigned int rightConnection = 0, leftConnection = 0;
	unsigned int postmp = 0;
	double disttmp = 0;
	pos = 0;
	for( typename IndexContainerType::const_iterator it = indexes.begin(); it != indexes.end(); ++it, ++i )
	{
		double tmp = this->Distance( *it, point, vector );
		if( tmp > threshold )
		{
			if( leftConnection == i )
			{
				leftConnection = i + 1;
			}
			else
			{
				rightConnection = i + 1;
				if( tmp > dist )
				{
					dist = tmp;
					pos = i;
				}
			}
		}
		else
		{
			if( rightConnection == i )
			{
				postmp = pos;
				disttmp = tmp;
			}
		}
	}
	if( 10 * rightConnection > 9 * i )
	{
		pos = postmp;
		dist = disttmp;
	}
	return pos;
}

template< class TGraphTraits >
typename SkeletonGraph< TGraphTraits >
::NodeIdentifierType SkeletonGraph< TGraphTraits >
::TryToSplit( NodePointerType node, double threshold, NodeIdentifierType & intersectionLabel )
{
	unsigned int pos;
	if( this->CanBeSplit( node->LineIndexes, node->Point, node->Vector, threshold, pos ) )
	{
		NodePointerType nNode = this->Split( node, pos, ++intersectionLabel );
		this->TryToSplit( node, threshold, intersectionLabel );
		this->TryToSplit( nNode, threshold, intersectionLabel );
	}
}

template< class TGraphTraits >
typename SkeletonGraph< TGraphTraits >
::NodePointerType SkeletonGraph< TGraphTraits >
::Split( NodePointerType node, unsigned int pos, NodeIdentifierType inter )
{
	NodePointerType nNode = this->CreateNewNode();
	nNode->BackLabel = node->BackLabel;
	nNode->FrontLabel = inter;
	node->BackLabel = inter;
	typename IndexContainerType::iterator idxit = node->LineIndexes.begin();
	typename ScalarContainerType::iterator scait = node->VesselnessValues.begin();
	for( unsigned int i = 0; i < pos; ++i )
	{
		++idxit;
		++scait;
	}
	nNode->LineIndexes.splice( nNode->LineIndexes.begin(), node->LineIndexes, idxit, node->LineIndexes.end() );
	nNode->VesselnessValues.splice( nNode->VesselnessValues.begin(), node->VesselnessValues, scait, node->VesselnessValues.end() );
	this->InitializeNode( nNode );
	this->CreateNewEdge( nNode->Identifier, nNode->Identifier );
	this->InitializeEdge( this->CreateNewEdge( nNode->Identifier, node->Identifier ), *nNode, *node );
	this->InitializeEdge( this->CreateNewEdge( node->Identifier, nNode->Identifier ), *node, *nNode );
	nNode->OutgoingEdgesFront.push_back( nNode->OutgoingEdges[1] );
	nNode->IncomingEdgesFront.push_back( nNode->IncomingEdges[1] );
	for( unsigned int i = node->OutgoingEdgesFront.size() + 1; i < node->OutgoingEdges.size() - 1; ++i )
	{
		EdgePointerType edge = this->GetEdgePointer( node->OutgoingEdges[i] );
		edge->SourceIdentifier = nNode->Identifier;
		this->InitializeEdge( edge, *nNode, this->GetTargetNode( node->OutgoingEdges[i] ) );
		edge = this->GetEdgePointer( node->IncomingEdges[i] );
		edge->TargetIdentifier = nNode->Identifier;
		this->InitializeEdge( edge, this->GetSourceNode( node->IncomingEdges[i] ), *nNode );
		nNode->OutgoingEdges.push_back( node->OutgoingEdges[i] );
		nNode->IncomingEdges.push_back( node->IncomingEdges[i] );
		nNode->OutgoingEdgesBack.push_back( node->OutgoingEdges[i] );
		nNode->IncomingEdgesBack.push_back( node->IncomingEdges[i] );
	}
	for( unsigned int i = node->OutgoingEdgesFront.size() + 1; i < node->OutgoingEdges.size() - 1; ++i )
	{
		this->InitializeNodeInteraction( this->GetTargetNodePointer( node->OutgoingEdges[i] ) );
	}
	for( unsigned int i = node->OutgoingEdgesFront.size() + 1; i < node->OutgoingEdges.size() - 1; )
	{
		node->OutgoingEdges.erase( node->OutgoingEdges.begin() + i );
		node->IncomingEdges.erase( node->IncomingEdges.begin() + i );
	}
	for( unsigned int i = 0; i < node->OutgoingEdgesBack.size(); )
	{
		node->OutgoingEdgesBack.erase( node->OutgoingEdgesBack.begin() + i );
		node->IncomingEdgesBack.erase( node->IncomingEdgesBack.begin() + i );
	}
	node->OutgoingEdgesBack.push_back( nNode->IncomingEdges[1] );
	node->IncomingEdgesBack.push_back( nNode->OutgoingEdges[1] );
	this->InitializeNodeInteraction( nNode );
	return nNode;
}

template< class TGraphTraits >
void SkeletonGraph< TGraphTraits >
::UpdateNodesWeightWithoutInteraction()
{
	for( typename NodeContainerType::Iterator it = this->GetNodeContainer()->Begin(); it != this->GetNodeContainer()->End(); ++it )
	{
		it.Value().Weight = this->GetFiberLikelihood( it.Value() );
	}
}

template< class TGraphTraits >
typename SkeletonGraph< TGraphTraits >
::NodeWeightType SkeletonGraph< TGraphTraits >
::GetFiberLikelihood( const NodeType & node ) const
{
	return node.MeanVesselness * this->GeometricLikelihood( node );
}

template< class TGraphTraits >
double SkeletonGraph< TGraphTraits >
::GeometricLikelihood( const NodeType & node ) const
{
	return std::log( node.LineIndexes.size() ) / ( 1 + std::abs( node.AverageDistance - node.TheoricalAverageDistance ) );
}


} // namespace itk

#endif
