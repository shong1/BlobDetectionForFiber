#ifndef __itkAVLIntervalTreeNode_hxx
#define __itkAVLIntervalTreeNode_hxx

#include "itkAVLIntervalTreeNode.h"

namespace itk
{

template< class TValue >
AVLIntervalTreeNode< TValue >::AVLIntervalTreeNode() : m_Low( NULL ), m_Up( NULL )
{
}

template< class TValue >
void AVLIntervalTreeNode< TValue >
::Set( Self *other )
{
	TreeNode< ValueType >::m_Data = other->Get();
	this->SetUp( other->GetUpLeaf() );
	this->SetLow( other->GetLowLeaf() );
}

template< class TValue >
void AVLIntervalTreeNode< TValue >
::Set( const ValueType & element )
{
	TreeNode< ValueType >::m_Data = element;
	this->SetLow( this );
	this->SetUp( this );
}

template< class TValue >
template< class TSelf >
TSelf* AVLIntervalTreeNode< TValue >
::GetRight( const ValueType element )
{
	return this->template GetNode< TSelf >( element )->template GetRight< TSelf >();
}

template< class TValue >
template< class TSelf >
TSelf* AVLIntervalTreeNode< TValue >
::GetRight()
{
	std::cout << " GetRight " << std::endl;
	if( this == NULL )
	{
		return NULL;
	}
	TSelf *node = dynamic_cast< TSelf* >( this );
	TSelf *tmp;
	unsigned int count = 0;
	std::cout << node->GetLow() << node->GetUp() << std::endl;
	do
	{
		tmp = node;
		node = node->GetParent();
		++count;
	}while( node != NULL && node->GetParent() != NULL && node->ChildPosition( tmp ) == Right );
	if( node != NULL && ( count == 1 || node->ChildPosition( tmp ) == Left ) )
	{
		return node->template GetChild< TSelf >( Right )->template GetLowLeaf< TSelf >();
	}
	else
	{
		return dynamic_cast< TSelf* >( this );
	}
}

template< class TValue >
template< class TSelf >
TSelf* AVLIntervalTreeNode< TValue >
::GetLeft( const ValueType element )
{
	return this->template GetNode< TSelf >( element )->template GetLeft< TSelf >();
}

template< class TValue >
template< class TSelf >
TSelf* AVLIntervalTreeNode< TValue >
::GetLeft()
{
	std::cout << " GetLeft " << std::endl;
	if( this == NULL )
	{
		return NULL;
	}
	TSelf *node = dynamic_cast< TSelf* >( this );
	TSelf *tmp;
	unsigned int count = 0;
	do
	{
		std::cout << " node " << node->GetLow() << node->GetUp() << std::endl;
		++count;
		tmp = node;
		node = node->GetParent();
	} while( node != NULL && node->GetParent() != NULL && node->ChildPosition( tmp ) == Left );
	std::cout << " hoho " << std::endl;
	if( node != NULL && ( count == 1 || node->ChildPosition( tmp ) == Right ) )
	{
		std::cout << " count " << count << " node " << node->GetLow() << node->GetUp() << std::endl;
		std::cout << node->GetParent() << std::endl;
		return node->template GetChild< TSelf >( Left )->template GetUpLeaf< TSelf >();
	}
	else
	{
		std::cout << " else GetLeft " << std::endl;
		return dynamic_cast< TSelf* >( this );
	}
}

template< class TValue >
template< class TSelf >
TSelf* AVLIntervalTreeNode< TValue >
::Add( TSelf *tmp )
{
	if( tmp == NULL )
	{
		return dynamic_cast< TSelf* >( this );
	}
	if( tmp->GetParent() != NULL )
	{
		tmp->GetParent()->AddChild( tmp->GetParent()->ChildPosition( tmp ), NULL );
		tmp->SetParent( NULL );
	}
	if( this == NULL )
	{
		return tmp;
	}
	tmp->SetUp( tmp );
	TSelf *node = dynamic_cast< TSelf* >( this ), *newRoot = dynamic_cast< TSelf* >( this );
	while( 1 )
	{
		if( *node <= *tmp )
		{
			if( !node->HasChildren() )
			{
				SmartPointer< TSelf > link = TSelf::New();
				link->Register();
				link->SetLow( node );
				link->SetUp( tmp );
				link->ShareInfo( dynamic_cast< TSelf* >( this ) );
				if( node->GetParent() != NULL )
				{
					node->GetParent()->AddChild( node->GetParent()->ChildPosition( node ), link );
				}
				else
				{
					newRoot = link;
				}
				link->AddChild( Left, node );
				link->AddChild( Right, tmp );
				node->UpdateDepthToRoot( this );
				newRoot =  newRoot->template Equilibrate< TSelf >( link );
				if( newRoot != NULL )
				return newRoot;
			}
			if( node->GetChild( Right ) == NULL )
			{
				node->SetUp( tmp );
				node->AddChild( Right, tmp );
				node->UpdateDepthToRoot( this );
				return newRoot->template Equilibrate< TSelf >( tmp );
			}
			//if( node->GetUp() < tmp->Get() )
			if( node->IsOutsideUpper( tmp ) )
			{
				node->SetUp( tmp );
			}
			node = node->GetChild( Right );
		}
		else
		{
			if( !node->HasChildren() )
			{
				SmartPointer< TSelf > link = TSelf::New();
				link->Register();
				link->SetLow( tmp );
				link->SetUp( node );
				link->ShareInfo( dynamic_cast< TSelf* >( this ) );
				if( node->GetParent() != NULL )
				{
					node->GetParent()->AddChild( node->GetParent()->ChildPosition( node ), link );
				}
				else
				{
					newRoot = link;
				}
				link->AddChild( Left, tmp );
				link->AddChild( Right, node );
				node->UpdateDepthToRoot( this );
				return newRoot->template Equilibrate< TSelf >( link );
			}
			if( node->GetChild( Left ) == NULL )
			{
				node->SetLow( tmp );
				node->AddChild( Left, tmp );
				node->UpdateDepthToRoot( this );
				return newRoot->template Equilibrate< TSelf >( tmp );
			}
			//if( node->GetLow() > tmp->Get() )
			if( node->IsOutsideLower( tmp ) )
			{
				node->SetLow( tmp );
			}
			node = node->template GetChild< TSelf >( Left );
		}
	}
}

template< class TValue >
template< class TSelf >
bool AVLIntervalTreeNode< TValue >
::IsOutsideUpper( TSelf* node ) const
{
	return this->GetUp() < node->GetLow();
}

template< class TValue >
template< class TSelf >
bool AVLIntervalTreeNode< TValue >
::IsOutsideLower( TSelf* node ) const
{
	return this->GetLow() > node->GetUp();
}

template< class TValue >
template< class TSelf >
TSelf *AVLIntervalTreeNode< TValue >
::Remove( TSelf *X )
{
	std::cout << " here " << std::endl;
	if( X == NULL )
	{
		return dynamic_cast< TSelf* >( this );
	}
	if( this == NULL )
	{
		return dynamic_cast< TSelf* >( this );
	}
	bool isGlobalMaximum = X->GetUp() == this->GetUp();
	bool isGlobalMinimum = X->GetLow() == this->GetLow();
	if( X->GetParent() != NULL )
	{
		//X->UnRegister();
	}
	TSelf *p = X->template GetParent< TSelf >(), *pp = NULL;
	TSelf *newRoot = dynamic_cast< TSelf* >( this );
	TSelf *leaf = NULL;
	if( p != NULL )
	{
		std::cout << " p != NULL " << std::endl;
		ChildIdentifier side = p->ChildPosition( X );
		std::cout << " ici " << std::endl;
		p->AddChild( side, NULL );
		std::cout << " ici " << std::endl;
		leaf = p->GetChild( !side );
		pp = p->template GetParent< TSelf >();
		if( pp != NULL )
		{
			std::cout << " but " << std::endl;
			pp->AddChild( pp->ChildPosition( p ), p->GetChild( !side ) );
		}
		else
		{
			std::cout << " else " << std::endl;
			newRoot = p->template GetChild< TSelf >( !side );
			newRoot->SetParent( NULL );
		}
	}
	else
	{
		return NULL;
	}
	std::cout << " here " << std::endl;
	if( isGlobalMaximum && pp != NULL )
	{
		ValueType up = leaf->GetUp();
		TSelf* node = pp;
		while( node != NULL )
		{
			node->SetUp( leaf );
			node = node->template GetParent< TSelf >();
		}
	}
	if( isGlobalMinimum && pp != NULL )
	{
		ValueType low = leaf->GetLow();
		TSelf* node = pp;
		while( node != NULL )
		{
			node->SetLow( leaf );
			node = node->template GetParent< TSelf >();
		}
	}
	X->SetParent( NULL );
	return newRoot->template Equilibrate< TSelf >( pp );
}

template< class TValue >
template< class TSelf >
TSelf* AVLIntervalTreeNode< TValue >
::Merge( TSelf *root )
{
	if( root == NULL )
	{
		std::cout << " root is NULL " << std::endl;
		return dynamic_cast< TSelf* >( this );
	}
	std::cout << " root is NOT NULL " << std::endl;
	if( root->GetParent() != NULL )
	{
		root->GetParent()->AddChild( root->GetParent()->ChildPosition( root ), NULL );
		root->SetParent( NULL );
	}
	if( this == NULL )
	{
		return root;
	}
	TSelf *newRoot = dynamic_cast< TSelf* >( this );
	if( root->GetDepth() > this->GetDepth() )
	{
		TSelf *tmp = dynamic_cast< TSelf* >( this );
		newRoot = root;
		root = tmp;
	}
	std::cout << " Merge newRoot " << *newRoot << std::endl;
	std::cout << " Merge root " << *root << std::endl;
	ChildIdentifier side;
	bool reverse = root->IsOutsideUpper( newRoot );
	if( reverse )
	{
		side = Left;
	}
	else
	{
		side = Right;
	}
	TSelf *r = newRoot;
	while( r->GetDepth() > root->GetDepth() )
	{
		if( side == Left )
		{
			r->SetLow( root->GetLowLeaf() );
		}
		else
		{
			r->SetUp( root->GetUpLeaf() );
		}
		r = r->GetChild( side );
	}
	SmartPointer< TSelf > link = TSelf::New();
	if( side == Right )
	{
		link->SetUp( root->GetUpLeaf() );
		link->SetLow( r->GetLowLeaf() );
	}
	else
	{
		link->SetLow( root->GetLowLeaf() );
		link->SetUp( r->GetUpLeaf() );
	}
	link->ShareInfo( dynamic_cast< TSelf* >( this ) );
	link->Register();
	if( r->GetParent() == NULL )
	{
		newRoot = link;
	}
	TSelf *rp = r->GetParent();
	rp->AddChild( side, link );
	link->AddChild( !side, r );
	link->AddChild( side, root );
	r = link;
	r->UpdateDepthToRoot( newRoot );		
	return newRoot->template Equilibrate< TSelf >( root, true );
}

template< class TValue >
template< class TSelf >
TSelf* AVLIntervalTreeNode< TValue >
::Split( TSelf *node, TSelf **root )
{
	if( node == NULL )
	{
		return dynamic_cast< TSelf* >( this );
	}
	ValueType element = node->Get();
	TSelf *tmp, *newRoot = NULL;
	node = dynamic_cast< TSelf* >( this );
	*root = NULL;
	std::cout << " split element " << element << std::endl;
	std::cout << " IICII " << std::endl;
	while( node != NULL )
	{
		std::cout << *node << std::endl;
		if( *node < element )// && !( node->GetChild( Left ) != NULL && *(node->GetChild( Left )) == element && node->GetChild( Right ) != NULL && *(node->GetChild( Right )) != element ) )
		{
			std::cout << "  Split Right " << std::endl;
			tmp = node->GetChild( Right );
			if( tmp != NULL )
			if( node->GetChild( Left ) == NULL )
			{
				*root = (*root)->Add( node );
				break;
			}
			else
			{
				*root = (*root)->Merge( node->GetChild( Left ) );
			}
		}
		else
		{
			std::cout << "  Split Left " << std::endl;
			tmp = node->GetChild( Left );
			if( node->GetChild( Right ) == NULL )
			{
				return newRoot->Add( node );
			}
			else
			{
				newRoot = newRoot->Merge( node->GetChild( Right ) );
			}
	if( *root != NULL )
	std::cout << " split *root " << (*root)->GetLow() << (*root)->GetUp() << std::endl;
	if( newRoot != NULL )
	std::cout << " split newRoot " << newRoot->GetLow() << newRoot->GetUp() << std::endl;
		}
		node->WipeChildren();
		if( node->HasChildren() )
		{
			node->UnRegister();
		}
		node = tmp;
		tmp->SetParent( NULL );
	}
	if( *root != NULL )
	std::cout << " split *root " << (*root)->GetLow() << (*root)->GetUp() << std::endl;
	if( newRoot != NULL )
	std::cout << " split newRoot " << newRoot->GetLow() << newRoot->GetUp() << std::endl;
	return newRoot;
}

template< class TValue >
template< class TSelf >
TSelf* AVLIntervalTreeNode< TValue >
::RotateRight( TSelf*  node )
{
	TSelf* newRoot = Superclass::template RotateRight< TSelf >( node );
	if( newRoot != NULL && node != NULL )
	{
		if( node->GetChild( Left ) != NULL )
		{
			node->SetLow( node->GetChild( Left )->GetLowLeaf() );
		}
		if( node->GetParent() != NULL )
		{
			node->GetParent()->SetUp( node->GetUpLeaf() );
		}
	}
	return newRoot;
}

template< class TValue >
template< class TSelf >
TSelf* AVLIntervalTreeNode< TValue >
::RotateLeft( TSelf* node )
{
	TSelf* newRoot = Superclass::template RotateLeft< TSelf >( node );
	if( newRoot != NULL && node != NULL )
	{
		if( node->GetChild( Right ) != NULL )
		{
			node->SetUp( node->GetChild( Right )->GetUpLeaf() );
		}
		if( node->GetParent() != NULL )
		{
			node->GetParent()->SetLow( node->GetLowLeaf() );
		}
	}
	return newRoot;
}

template< class TValue >
template< class TSelf >
TSelf* AVLIntervalTreeNode< TValue >
::Equilibrate( TSelf *node, bool recursive )
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
			newRoot = newRoot->template RotateLeft< TSelf >( L );
			newRoot = newRoot->template RotateRight< TSelf >( node );
		}
		if( !recursive )
		{
			return newRoot;
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
			newRoot = newRoot->template RotateRight< TSelf >( R );
			newRoot = newRoot->template RotateLeft< TSelf >( node );
		}
		if( !recursive )
		{
			return newRoot;
		}
	}
	return newRoot->template Equilibrate< TSelf >( node->GetParent(), recursive );
}

template< class TValue >
template< class TSelf >
bool AVLIntervalTreeNode< TValue >
::operator<( const TSelf &other )
{
	if( this->GetChild( Left ) != NULL )
	{
		return this->GetChild( Left )->GetUp() <= other.GetUp();
	}
	if( this->GetChild( Right ) != NULL )
	{
		return this->GetChild( Right )->GetLow() >= other.GetLow();
	}
	return this->GetUp() < other.GetLow();
}

template< class TValue >
template< class TSelf >
bool AVLIntervalTreeNode< TValue >
::operator<=( const TSelf &other )
{
	if( this->GetChild( Left ) != NULL )
	{
		return this->GetChild( Left )->GetUp() <= other.GetUp();
	}
	if( this->GetChild( Right ) != NULL )
	{
		return this->GetChild( Right )->GetLow() >= other.GetLow();
	}
	return this->GetUp() <= other.GetLow();
}

template< class TValue >
template< class TSelf >
bool AVLIntervalTreeNode< TValue >
::operator>( const TSelf &other )
{
	if( this->GetChild( Right ) != NULL )
	{
		return this->GetChild( Right )->GetLow() <= other.GetLow();
	}
	if( this->GetChild( Left ) != NULL )
	{
		return this->GetChild( Left )->GetUp() > other.GetUp();
	}
	return this->GetUp() > other.GetLow();
}

template< class TValue >
template< class TSelf >
bool AVLIntervalTreeNode< TValue >
::operator>=( const TSelf &other )
{
	if( this->GetChild( Right ) != NULL )
	{
		return this->GetChild( Right )->GetLow() <= other.GetLow();
	}
	if( this->GetChild( Left ) != NULL )
	{
		return this->GetChild( Left )->GetUp() > other.GetUp();
	}
	return this->GetUp() >= other.GetLow();
}

template< class TValue >
template< class TSelf >
bool AVLIntervalTreeNode< TValue >
::operator==( const TSelf &other )
{
	return this->GetUp() == other.GetUp() && this->GetLow() == other.GetLow();
}

template< class TValue >
template< class TSelf >
bool AVLIntervalTreeNode< TValue >
::operator!=( const TSelf &other )
{
	return this->GetUp() != other.GetUp() || this->GetLow() != other.GetLow();
}

template< class TValue >
template< class TSelf >
bool AVLIntervalTreeNode< TValue >
::operator<( const ValueType & element )
{
	if( this->GetChild( Left ) != NULL )
	{
		return this->GetChild( Left )->GetUp() <= element;
	}
	if( this->GetChild( Right ) != NULL )
	{
		return this->GetChild( Right )->GetLow() > element;
	}
	return this->GetUp() < element;
}

template< class TValue >
template< class TSelf >
bool AVLIntervalTreeNode< TValue >
::operator<=( const ValueType & element )
{
	if( this->GetChild( Left ) != NULL )
	{
		return this->GetChild( Left )->GetUp() <= element;
	}
	if( this->GetChild( Right ) != NULL )
	{
		return this->GetChild( Right )->GetLow() > element;
	}
	return this->GetUp() <= element;
}

template< class TValue >
template< class TSelf >
bool AVLIntervalTreeNode< TValue >
::operator>( const ValueType & element )
{
	if( this->GetChild( Right ) != NULL )
	{
		return this->GetChild( Right )->GetLow() <= element;
	}
	if( this->GetChild( Left ) != NULL )
	{
		return this->GetChild( Left )->GetUp() > element;
	}
	return this->GetUp() > element;
}

template< class TValue >
template< class TSelf >
bool AVLIntervalTreeNode< TValue >
::operator>=( const ValueType & element )
{
	if( this->GetChild( Right ) != NULL )
	{
		return this->GetChild( Right )->GetLow() <= element;
	}
	if( this->GetChild( Left ) != NULL )
	{
		return this->GetChild( Left )->GetUp() > element;
	}
	return this->GetUp() >= element;
}

template< class TValue >
template< class TSelf >
bool AVLIntervalTreeNode< TValue >
::operator==( const ValueType & element )
{
	return this->GetUp() == element && this->GetLow() == element;
}

template< class TValue >
template< class TSelf >
bool AVLIntervalTreeNode< TValue >
::operator!=( const ValueType & element )
{
	return this->GetUp() != element || this->GetLow() != element;
}

template< class T >
std::ostream& operator<<( std::ostream& flux, AVLIntervalTreeNode< T > const& node )
{
	flux << node.GetLow() << " " << node.GetUp() << std::endl;
	if( node.GetParent() != NULL )
	{
		flux << " parent " << node.GetParent()->GetLow() << " " << node.GetParent()->GetUp() << std::endl;
	}
	if( node.GetChild( 0 ) != NULL )
	{
		flux << " Left " << std::endl;
		flux << *(node.GetChild( 0 ));
		if( node.GetParent() != NULL )
		{
			flux << " parent " << node.GetParent()->GetLow() << " " << node.GetParent()->GetUp() << std::endl;
		}
	}
	if( node.GetChild( 1 ) != NULL )
	{
		flux << " Right " << std::endl;
		flux << *(node.GetChild( 1 ));
	}
	return flux;
}

}// namespace itk

#endif
