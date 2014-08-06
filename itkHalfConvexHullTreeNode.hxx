#ifndef __itkHalfConvexHullTreeNode_hxx
#define __itkHalfConvexHullTreeNode_hxx

#include "itkHalfConvexHullTreeNode.h"

namespace itk
{

template< class TVector >
bool HalfConvexHullTreeNode< TVector >
::ChildrenAreEmpty() const
{
	return this->GetDepth() == 2 && this->GetChild( Left )->GetQRoot() == NULL && this->GetChild( Right )->GetQRoot() == NULL;
}

template< class TVector >
template< class TSelf >
TSelf* HalfConvexHullTreeNode< TVector >
::Add( ConcatenableQueueNodeType * node, TSelf* newRoot )
{
	TSelf *current = dynamic_cast< TSelf* >( this );
	if( newRoot == NULL )
	{
		newRoot = current;
	}
	// Initializing case
	if( this == NULL && node != NULL )
	{
		SmartPointer< TSelf > p = TSelf::New();
		current = p.GetPointer();
		current->Register();
		SmartPointer< ConcatenableQueueType > q = ConcatenableQueueType::New();
		q.GetPointer()->Register();
		current->SetQ( q );
		current->GetQ()->Add( node );
		return current;
	}
	else if( this == NULL && node == NULL )
	{
		return NULL;
	}
	else if( this->GetQRoot() == NULL && node != NULL )
	{
		current->GetQ()->SetRoot( node );
		current->SetBridge( node );
		return newRoot->Up( current );
	}
	bool nodeAdded = false;

	while( 1 )
	{
		// If there are no children, we have to create them
		// because we go deeper in the tree at each iteration
		if( !current->HasChildren() )
		{
			SmartPointer< ConcatenableQueueType > nq = ConcatenableQueueType::New();
			ConcatenableQueueType *newQueue = nq.GetPointer();
			newQueue->Register();
			SmartPointer< ConcatenableQueueType > onq = ConcatenableQueueType::New();
			ConcatenableQueueType *otherNewQueue = onq.GetPointer();
			otherNewQueue->Register();
			typename TSelf::Pointer r = TSelf::New();
			TSelf *right = r.GetPointer();
			right->Register();
			right->SetTree( this->GetTree() );
			right->SetQ( newQueue );
			current->AddChild( Right, right );
			right->UpdateDepthToRoot( this );
			typename TSelf::Pointer l = TSelf::New();
			TSelf *left = l.GetPointer();
			left->Register();
			left->SetTree( this->GetTree() );
			left->SetQ( otherNewQueue );
			current->AddChild( Left, left );
			if( current->GetQRoot() != NULL )
			{
				// If we already have a node with same element in Dimension,
				// we store it into an array at the location of the node.
				// If not, it messes the algorithm.
				ConcatenableQueueNodeType *dup = current->GetQRoot()->GetNode( node->Get() );
				if( dup != NULL )
				{
					if( dup != dup->Duplicate( node, m_Tree->GetContour() ) )
					{
						dup = node;
						if( dup->GetParent() == NULL )
						{
							current->GetQ()->SetRoot( dup );
						}
						if( current->GetQRoot()->GetSize() == 1 )
						{
							if( current->GetParent() != NULL )
							{
								return newRoot->Up( current->GetParent() );
							}
							return current;
						}
						current->GetQ()->Split( node, current->GetChild( Left )->GetQ() );
						current->GetChild( Right )->GetQ()->SetRoot( current->GetQRoot() );
						current->GetQ()->SetRoot( NULL );
						return newRoot->Up( current );
					}
					if( current->GetParent() != NULL )
					{
						return newRoot->Up( current->GetParent() );
					}
					return current;
				}
				// We basically have to add a node in this case
				// Because no children (se above) means we cannot go deeper
				if( node->GetLow()[m_Tree->GetDimension()] > current->GetQRoot()->GetUp()[m_Tree->GetDimension()]
						|| node->GetLow()[m_Tree->GetDimension()] == current->GetQRoot()->GetUp()[m_Tree->GetDimension()]
						&& node->GetLow()[m_Tree->GetOtherDimension()] >= current->GetQRoot()->GetUp()[m_Tree->GetOtherDimension()] )
				{
					current->GetChild( Right )->GetQ()->Add( node );
					current->GetChild( Left )->GetQ()->SetRoot( current->GetQRoot() );
					current->GetQ()->SetRoot( NULL );
					return newRoot->Up( current );
				}
				else if( node->GetUp()[m_Tree->GetDimension()] <= current->GetQRoot()->GetLow()[m_Tree->GetDimension()]
						|| node->GetUp()[m_Tree->GetDimension()] == current->GetQRoot()->GetLow()[m_Tree->GetDimension()]
						&& node->GetUp()[m_Tree->GetOtherDimension()] <= current->GetQRoot()->GetLow()[m_Tree->GetOtherDimension()] )
				{
					current->GetChild( Left )->GetQ()->Add( node );
					current->GetChild( Right )->GetQ()->SetRoot( current->GetQRoot() );
					current->GetQ()->SetRoot( NULL );
					return newRoot->Up( current );
				}
			}
			// Empty queue : we create it
			else
			{
				current->GetQ()->Add( node );
				return newRoot->Up( current );
			}
		}
		// Children exist but have no queue in it
		if( current->ChildrenAreEmpty() )
		{
			// If we are outside upper : it means upper than any points of the set
			// the node is in the hull
			if( current->GetQRoot()->IsOutsideUpper( node ) )
			{
				current->GetChild( Right )->GetQ()->Add( node );
				current->GetChild( Left )->GetQ()->SetRoot( current->GetQRoot() );
				current->GetQ()->SetRoot( NULL );
				return newRoot->Up( current );
			}
			// If we are outside lower : it means lower than any points of the set
			// the node is in the hull
			else if( current->GetQRoot()->IsOutsideLower( node ) )
			{
				current->GetChild( Left )->GetQ()->Add( node );
				current->GetChild( Right )->GetQ()->SetRoot( current->GetQRoot() );
				current->GetQ()->SetRoot( NULL );
				return newRoot->Up( current );
			}
			// else, we have to go down
			else
			{
				// checking if there is a duplicate
				ConcatenableQueueNodeType *dup = current->GetQRoot()->GetNode( node->Get() );
				if( dup != NULL )
				{
					if( dup != dup->Duplicate( node, m_Tree->GetContour() ) )
					{
						dup = node;
						if( dup->GetParent() == NULL )
						{
							current->GetQ()->SetRoot( dup );
						}
						if( current->GetQRoot()->GetSize() == 1 )
						{
							if( current->GetParent() != NULL )
							{
								return this->Up( current->GetParent() );
							}
							return current;
						}
						current->GetQ()->Split( node, current->GetChild( Left )->GetQ() );
						current->GetChild( Right )->GetQ()->SetRoot( current->GetQRoot() );
						current->GetQ()->SetRoot( NULL );
						return newRoot->Up( current );
					}
					if( current->GetParent() != NULL )
					{
						return newRoot->Up( current->GetParent() );
					}
					return current;
				}
				// We split at the location of the new node, and we go down
				current->GetQ()->Split( node, current->GetChild( Left )->GetQ() );
				current->GetChild( Right )->GetQ()->SetRoot( current->GetQRoot() );
				current->GetQ()->SetRoot( NULL );
				current = current->GetChild( Left );
			}
		}
		// The children are not empty. We decide which way to go, and we go down
		else
		{
			current->Down( current->GetChild( Right ) );
			if( current->GetChild( Right )->GetQRoot() == NULL 
					&& current->GetChild( Left )->GetQRoot() != NULL 
					&& current->GetChild( Left )->GetQRoot()->GetUp()[m_Tree->GetDimension()] < node->Get()[m_Tree->GetDimension()]
					|| current->GetChild( Right )->GetQRoot() != NULL
					&& current->GetChild( Right )->GetQRoot()->GetLow()[m_Tree->GetDimension()] <= node->Get()[m_Tree->GetDimension()] )
			{
				current = current->GetChild( Right );
			}
			else
			{
				current = current->GetChild( Left );
			}
		}
	}
}

template< class TVector >
template< class TSelf >
TSelf* HalfConvexHullTreeNode< TVector >
::Remove( const VectorType & element )
{
	TSelf* current = dynamic_cast< TSelf* >( this );
	TSelf* newRoot = current;
	while( 1 )
	{
		// Children? we go down (we have to get to the leaves)
		if( current->HasChildren() && current->GetBridge() != NULL )
		{
			current->Down( current->GetChild( Right ) );
			if( current->GetChild( Right )->GetQRoot()->GetLow()[m_Tree->GetDimension()] > element[m_Tree->GetDimension()] )
			{
				current = current->GetChild( Left );
			}
			else
			{
				current = current->GetChild( Right );
			}
		}
		// No children? time to remove
		else
		{
			ConcatenableQueueNodeType *node = current->GetQRoot()->GetNode( element );
			if( node != NULL )
			{
				if( node->Get() == element )
				{
					// No duplicate : easy
					if( node->GetDuplicate() == NULL )
					{
						current->GetQ()->Remove( node );
						if( current->GetQRoot() == NULL )
						{
							TSelf* p = current->GetParent();
							if( p != NULL )
							{
								ChildIdentifier pos = p->ChildPosition( current );
								p->GetChild( !pos )->SetParent( p->GetParent() );
								TSelf* pp = p->GetParent();
								if( pp != NULL )
								{
									pp->AddChild( pp->ChildPosition( p ), p->GetChild( !pos ) ); 
									pp->UpdateDepthToRoot( newRoot );
									newRoot = newRoot->Equilibrate( pp );
									return newRoot->Up( pp );
								}
								return newRoot->Up( p );
							}
							return newRoot;
						}
					}
					// Duplicate : remove and update it
					else
					{
						ConcatenableQueueNodeType* tmp = node, *buf = tmp;;
						while( tmp != NULL && tmp->Get() != element )
						{
							buf = tmp;
							tmp = tmp->GetDuplicate();
						}
						if( tmp == NULL )
						{
							return newRoot->Up( current );
						}
						if( tmp != node )
						{
							buf->SetDuplicate( tmp->GetDuplicate() );
							return newRoot->Up( current );
						}
						ConcatenableQueueNodeType* dup = node->GetDuplicate();
						dup->SetParent( NULL );
						current->GetQ()->Remove( node );
						return current->Add( dup, newRoot );
					}
					// If we erased the last point of the current queue
					if( current->GetQRoot() == NULL )
					{
						current->GetParent()->GetQ()->SetRoot( current->GetParent()->GetChild( !current->GetParent()->ChildPosition( current ) ) );
						current = current->GetParent();
						current->SetBridge( NULL );
					}
				}
			}
			return newRoot->Up( current );
		}
	}
}

template< class TVector >
template< class TSelf >
TSelf* HalfConvexHullTreeNode< TVector >
::Down( TSelf *to )
{
	// See paper
	if( this == NULL || to == NULL )
	{
		return dynamic_cast< TSelf* >( this );
	}
	TSelf *current = dynamic_cast< TSelf* >( this );
	ConcatenableQueueNodeType *interval = to->GetQRoot();
	while( 1 )
	{
		ChildIdentifier side;
		if( current == to )
		{
			return to;
		}
		if( current->GetQRoot()->GetSize() == 1 )
		{
			current->SetBridge( current->GetQRoot() );
		}
		if( current->HasChildren() )
		{
			TSelf * tmp = to;
			bool keepGoing = true;
			while( interval == NULL && keepGoing  )
			{
				if( tmp->GetParent() == current )
				{
					side = current->ChildPosition( tmp );
					keepGoing = false;
				}
				else
				{
					tmp = tmp->GetParent();
					interval = tmp->GetQRoot();
				}
			}
			if( interval != NULL )
			{
				if( current->GetChild( Left ) != NULL && current->GetChild( Left )->GetQRoot() != NULL && current->GetChild( Left )->GetQRoot()->GetUp()[m_Tree->GetDimension()] >= interval->Get()[m_Tree->GetDimension()] )
				{
					side = Left;
				}
				else if( current->GetChild( Right ) != NULL && current->GetChild( Right )->GetQRoot() != NULL && current->GetChild( Right )->GetQRoot()->GetLow()[m_Tree->GetDimension()] <= interval->Get()[m_Tree->GetDimension()] )
				{
					side = Right;
				}
				else
				{
					this->Up( current );
					return dynamic_cast< TSelf* >( this );
				}
			}
		}
		else
		{
			this->Up( current );
			return dynamic_cast< TSelf* >( this );
		}
		ConcatenableQueueNodeType* leftRoot = NULL;
		current->GetQ()->Split( current->GetBridge(), &leftRoot );
		current->GetChild( Right )->GetQ()->Merge( current->GetQ() );
		current->GetChild( Left )->GetQ()->Merge( leftRoot );
		current = current->GetChild( side );
		TSelf *tmp = current->GetParent();
	}
}

template< class TVector >
template< class TSelf >
TSelf* HalfConvexHullTreeNode< TVector >
::Up( TSelf *node, int number )
{
	// See paper
	if( this == NULL )
	{
		return NULL;
	}
	TSelf *newRoot = dynamic_cast< TSelf* >( this );
	while( node != NULL && number != 0 )
	{
		newRoot = newRoot->Equilibrate( node );

		TSelf *lson = node->GetChild( Left ), *rson = node->GetChild( Right );
		if( lson != NULL && rson != NULL )
		{
			ConcatenableQueueNodeType *bridge = NULL, *brjdge;
			this->UpdateBridge( lson->GetQRoot(), rson->GetQRoot(), &bridge, &brjdge ); 
			node->SetBridge( bridge );	
			ConcatenableQueueType *l = node->GetQ();
			ConcatenableQueueNodeType *root = NULL;
			if( bridge->GetRight() != bridge )
			{
				lson->GetQ()->Split( bridge, l );
			}
			else
			{
				node->GetQ()->Merge( lson->GetQ() );
			}
			ConcatenableQueueNodeType *leftNode = brjdge->GetLeft();
			if( leftNode != brjdge )
			{
				rson->GetQ()->Split( leftNode, &root );
			}
			l->Merge( rson->GetQ() );
			rson->GetQ()->SetRoot( root );
			ConcatenableQueueType *tmp = rson->GetQ();
		}
		node = node->GetParent();
		--number;
	}
	if( newRoot->GetQRoot() != NULL )
	return newRoot;
}

template< class TVector >
template< class TSelf >
TSelf* HalfConvexHullTreeNode< TVector >
::Equilibrate( TSelf *node )
{
	if( node == NULL )
	{
		return dynamic_cast< TSelf* >( this );
	}
	if( this == NULL )
	{
		return dynamic_cast< TSelf* >( this );
	}
	TSelf *L = node->template GetChild< TSelf >( Left ), *R = node->template GetChild< TSelf >( Right );
	TSelf *newRoot = dynamic_cast< TSelf* >( this );
	int balance = L->GetDepth() - R->GetDepth();
	if( balance == 2 )
	{
		Self *LL = L->template GetChild< TSelf >( Left ), *LR = L->template GetChild< TSelf >( Right );
		if( LL->GetDepth() > LR->GetDepth() )
		{
			newRoot = newRoot->template RotateRight< TSelf >( node );
		}
		else
		{
			L->Down( L->GetChild( Left ) );
			newRoot = newRoot->template RotateLeft< TSelf >( L );
			newRoot = newRoot->Up( L, 1 );
			newRoot = newRoot->template RotateRight< TSelf >( node );
		}
	}
	if( balance == -2 )
	{
		Self *RL = R->template GetChild< TSelf >( Left ), *RR = R->template GetChild< TSelf >( Right );
		if( RR->GetDepth() > RL->GetDepth() )
		{
			newRoot = newRoot->template RotateLeft< TSelf >( node );
		}
		else
		{
			R->Down( R->GetChild( Left ) );
			newRoot = newRoot->template RotateRight< TSelf >( R );
			newRoot = newRoot->Up( R, 1 );
			newRoot = newRoot->template RotateLeft< TSelf >( node );
		}
	}
	return newRoot;
}

template< class TVector >
template< class TSelf >
TSelf* HalfConvexHullTreeNode< TVector >
::RotateRight( TSelf *node )
{
	if( node == NULL )
	{
		return dynamic_cast< TSelf* >( this );
	}
	node->GetChild( Left )->Down( node->GetChild( Left )->GetChild( Left ) );
	return Superclass::template RotateRight< TSelf >( node );
}

template< class TVector >
template< class TSelf >
TSelf* HalfConvexHullTreeNode< TVector >
::RotateLeft( TSelf *node )
{
	if( node == NULL )
	{
		return dynamic_cast< TSelf* >( this );
	}
	if( node == NULL )
	{
		return dynamic_cast< TSelf* >( this );
	}
	node->GetChild( Right )->Down( node->GetChild( Right )->GetChild( Right ) );
	return  Superclass::template RotateLeft< TSelf >( node );
}

template< class TVector >
bool HalfConvexHullTreeNode< TVector >
::ContourSide( const VectorType & down, const VectorType & up, const VectorType & x, ScalarType * result ) const
{
	VectorType vec = up - down;
	ScalarType c =  vec[m_Tree->GetDimension()] * up[m_Tree->GetOtherDimension()] - vec[m_Tree->GetOtherDimension()] * up[m_Tree->GetDimension()];
	*result = x[m_Tree->GetDimension()] * vec[m_Tree->GetOtherDimension()] - x[m_Tree->GetOtherDimension()] * vec[m_Tree->GetDimension()] + c ;
	if( *result == 0 && x[m_Tree->GetDimension()] < up[m_Tree->GetDimension()] && x[m_Tree->GetDimension()] > down[m_Tree->GetDimension()] )
	{
		return m_Tree->GetContour();
	}
	if( *result > 0 && m_Tree->GetContour() == Left || *result < 0 && m_Tree->GetContour() == Right )
	{
		return m_Tree->GetContour();
	}
	return !m_Tree->GetContour();
}

// i < j
// ref : www.ics.uci.edu/~goodrich/teach/cs164/notes/sdarticle-37.pdf
template< class TVector >
void HalfConvexHullTreeNode< TVector >
::UpdateBridge( ConcatenableQueueNodeType *i, ConcatenableQueueNodeType *j, ConcatenableQueueNodeType **bridge, ConcatenableQueueNodeType **brjdge )
{
	if( i == NULL )
	{
		*bridge = NULL;
	}
	if( j == NULL )
	{
		*brjdge = NULL;
	}
	if( i == NULL && j != NULL )
	{
		*brjdge = j->GetLowLeaf();
	}
	if( i != NULL && j == NULL )
	{
		*bridge = i->GetUpLeaf();
	}
	if( i == NULL || j == NULL )
	{
		return;
	}
	// left i, right i, left j, right j, current j, current i
	ConcatenableQueueNodeType *li, *ri, *lj, *rj, *tj, *ti;

	std::cout << i << " " << j << std::endl;

	if( i->GetLow() != i->GetUp() )
	{
		li = i->GetChild( Left )->GetUpLeaf();
		ri = i->GetChild( Right )->GetLowLeaf();
	}
	else
	{
		li = i;
		ri = i;
	}
	if( j->GetLow() != j->GetUp() )
	{
		lj = j->GetChild( Left )->GetUpLeaf();
		rj = j->GetChild( Right )->GetLowLeaf();
	}
	else
	{
		lj = j;
		rj = j;
	}
	ti = li;
	tj = lj;

	// To not be stuck anywere in the routine
	// It is used so once we do i, once we do j
	bool state = true;
	ScalarType result;

	while( 1 )
	{
		ConcatenableQueueNodeType *lti, *rti, *ltj, *rtj;

		bool blti = !m_Tree->GetContour(), brti = m_Tree->GetContour(), bltj = !m_Tree->GetContour(), brtj = m_Tree->GetContour();

		// left and right of the current i and j, and their side
		lti = ti->GetLeft();
		rti = ti->GetRight();
		blti = this->ContourSide( ti->Get(), tj->Get(), lti->Get(), &result );
		brti = this->ContourSide( ti->Get(), tj->Get(), rti->Get(), &result );

		ltj = tj->GetLeft();
		rtj = tj->GetRight();
		bltj = this->ContourSide( ti->Get(), tj->Get(), ltj->Get(), &result );
		brtj = this->ContourSide( ti->Get(), tj->Get(), rtj->Get(), &result );

		// case a
		if( blti == brti && brti == bltj && bltj == brtj && blti != m_Tree->GetContour() )
		{
			*bridge = ti;
			*brjdge = tj;
			return;
		}

		// case b
		else if( bltj == brtj && blti == !bltj && brti == bltj )
		{
			if( li != ri )
			{
				if( ti == li )
				{
					i = i->GetChild( Left );
					if( i->HasChildren() )
					{
						li = i->GetChild( Left )->GetUpLeaf();
						ri = i->GetChild( Right )->GetLowLeaf();
					}
					else
					{
						li = i;
						ri = i;
					}
					ti = li;
				}
				else
				{
					if( tj == lj && state )
					{
						ti = li;
					}
				}
			}
			if( lj != rj )
			{
				if( tj == lj )
				{
					if( ti == ri )
					{
						tj = rj;
					}
				}
				else
				{
					j = j->GetChild( Right );
					if( j->HasChildren() )
					{
						lj = j->GetChild( Left )->GetUpLeaf();
						rj = j->GetChild( Right )->GetLowLeaf();
					}
					else
					{
						lj = j;
						rj = j;
					}
					tj = lj;
				}
			}
		}

		// case c
		else if( bltj == brtj && blti == bltj && brti == !bltj )
		{
			if( li != ri )
			{
				if( ti == li )
				{
					if( tj == lj && state )
					{
						ti = ri;
					}
				}
				else
				{
					i = i->GetChild( Right );
					if( i->HasChildren() )
					{
						li = i->GetChild( Left )->GetUpLeaf();
						ri = i->GetChild( Right )->GetLowLeaf();
					}
					else
					{
						li = i;
						ri = i;
					}
					ti = li;
				}
			}
			if( lj != rj )
			{
				if( tj == lj )
				{
					if( ti == li && state )
					{
						tj = rj;
					}
				}
				else
				{
					j = j->GetChild( Right );
					if( j->HasChildren() )
					{
						lj = j->GetChild( Left )->GetUpLeaf();
						rj = j->GetChild( Right )->GetLowLeaf();
					}
					else
					{
						lj = j;
						rj = j;
					}
					tj = lj;
				}
			}
		}

		// case d
		else if( blti == brti && brti == bltj && brti == !brtj )
		{
			if( li != ri )
			{
				if( ti == li )
				{
					i = i->GetChild( Left );
					if( i->HasChildren() )
					{
						li = i->GetChild( Left )->GetUpLeaf();
						ri = i->GetChild( Right )->GetLowLeaf();
					}
					else
					{
						li = i;
						ri = i;
					}
					ti = li;
				}
				else
				{
					if( tj == lj && state )
					{
						ti = li;
					}
				}
			}
			if( lj != rj )
			{
				if( tj == lj )
				{
					if( ti == ri )
					{
						tj = rj;
					}
				}
				else
				{
					j = j->GetChild( Right );
					if( j->HasChildren() )
					{
						lj = j->GetChild( Left )->GetUpLeaf();
						rj = j->GetChild( Right )->GetLowLeaf();
					}
					else
					{
						lj = j;
						rj = j;
					}
					tj = lj;
				}
			}
		}

		// case e
		else if( blti == brti && brti == !bltj && brti == brtj )
		{
			if( li != ri )
			{
				if( ti == li )
				{
					i = i->GetChild( Left );
					if( i->HasChildren() )
					{
						li = i->GetChild( Left )->GetUpLeaf();
						ri = i->GetChild( Right )->GetLowLeaf();
					}
					else
					{
						li = i;
						ri = i;
					}
					ti = li;
				}
				else
				{
					if( tj == rj && state )
					{
						ti = li;
					}
				}
			}
			if( lj != rj )
			{
				if( tj == lj )
				{
					j = j->GetChild( Left );
					if( j->HasChildren() )
					{
						lj = j->GetChild( Left )->GetUpLeaf();
						rj = j->GetChild( Right )->GetLowLeaf();
					}
					else
					{
						lj = j;
						rj = j;
					}
					tj = lj;
				}	
				else
				{
					if( ti == ri )
					{
						tj = lj;
					}
				}
			}
		}

		// case f
		else if( blti == m_Tree->GetContour() && bltj == brti && brtj == blti && brti != blti )
		{
			if( ri != li )
			{
				if( ti == li )
				{
					i = i->GetChild( Left );
					if( i->HasChildren() )
					{
						li = i->GetChild( Left )->GetUpLeaf();
						ri = i->GetChild( Right )->GetLowLeaf();
					}
					else
					{
						li = i;
						ri = i;
					}
					ti = li;
				}
				else
				{
					if( tj == lj && state )
					{
						ti = li;
					}
				}
			}
			if( rj != lj )
			{
				if( tj == lj )
				{
					if( ti == ri )
					{
						tj = rj;
					}
				}
				else
				{
					j = j->GetChild( Right );
					if( j->HasChildren() )
					{
						lj = j->GetChild( Left )->GetUpLeaf();
						rj = j->GetChild( Right )->GetLowLeaf();
					}
					else
					{
						lj = j;
						rj = j;
					}
					tj = lj;
				}
			}
		}

		// case g
		else if(  blti == m_Tree->GetContour() && brti == brtj && blti == bltj && blti != brti )
		{
			if( li != ri )
			{
				if( ti == li )
				{
					i = i->GetChild( Left );
					if( i->HasChildren() )
					{
						li = i->GetChild( Left )->GetUpLeaf();
						ri = i->GetChild( Right )->GetLowLeaf();
					}
					else
					{
						li = i;
						ri = i;
					}
					ti = li;
				}
				else
				{
					ti = li;
				}
			}
			state = !state;
		}

		// case h
		else if(  brti == m_Tree->GetContour() && brti == brtj && blti == bltj && blti != brti )
		{
			if( lj != rj )
			{
				if( tj == lj )
				{
					tj = rj;
				}
				else
				{
					j = j->GetChild( Right );
					if( j->HasChildren() )
					{
						lj = j->GetChild( Left )->GetUpLeaf();
						rj = j->GetChild( Right )->GetLowLeaf();
					}
					else
					{
						lj = j;
						rj = j;
					}
					tj = lj;
				}
			}
			state = !state;
		}

		// case i
		else if( brti == m_Tree->GetContour() && bltj == brti && brtj == blti && brti != blti )
		{
			VectorType tanti = rti->Get() - lti->Get();
			VectorType tantj= rtj->Get() - ltj->Get();
			float t = tantj[m_Tree->GetOtherDimension()] * ( ti->Get()[m_Tree->GetDimension()] - tj->Get()[m_Tree->GetDimension()] ) - tantj[m_Tree->GetDimension()] * ( ti->Get()[m_Tree->GetOtherDimension()] - tj->Get()[m_Tree->GetOtherDimension()] ) / (double) (tantj[m_Tree->GetOtherDimension()] * tanti[m_Tree->GetDimension()] - tantj[m_Tree->GetDimension()] * tanti[m_Tree->GetOtherDimension()] );

			// case i2
			if( !( ( ti->Get() + tanti * t )[m_Tree->GetDimension()] < j->GetLow()[m_Tree->GetDimension()] ) )
			{
				if( rj != lj )
				{
					if( tj == lj )
					{
						j = j->GetChild( Left );
						if( j->HasChildren() )
						{
							lj = j->GetChild( Left )->GetUpLeaf();
							rj = j->GetChild( Right )->GetLowLeaf();
						}
						else
						{
							lj = j;
							rj = j;
						}
						tj = lj;
					}
					else
					{
						tj = lj;
					}
				}
			}

			// case i1
			else
			{
				if( ri != li )
				{
					if( ti == li )
					{
						ti = ri;
					}
					else
					{
						i = i->GetChild( Right );
						if( i->HasChildren() )
						{
							li = i->GetChild( Left )->GetUpLeaf();
							ri = i->GetChild( Right )->GetLowLeaf();
						}
						else
						{
							li = i;
							ri = i;
						}
						ti = li;
					}
				}
			}
			state = !state;
		}
		
		// error
		else
		{
			if( ti == li && ti != ri )
			{
				std::cout << " error 1 " << std::endl;
				ti = ri;
			}
			else if( tj == lj && tj != rj )
			{
				std::cout << " error 2 " << std::endl;
				tj = rj;
			}
			else if( ti == ri && ti != li )
			{
				std::cout << " error 3 " << std::endl;
				ti = li;
			}
			else if( tj == rj && tj != lj)
			{
				std::cout << " error 4 " << std::endl;
				tj = lj;
			}
			else
			{
				std::cout << " error 5 " << std::endl;
			}
		}
		state = !state;
	}
}

}// namespace itk

#endif
