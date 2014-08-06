#ifndef __itkAVLTreeNode_hxx
#define __itkAVLTreeNode_hxx

#include "itkAVLTreeNode.h"

namespace itk
{

template< class TValue >
AVLTreeNode< TValue >::AVLTreeNode()
{
}

template< class TValue >
template< class TSelf >
TSelf* AVLTreeNode< TValue >
::Add( TSelf *tmp )
{
	TSelf *newRoot = Superclass::template Add< TSelf >( tmp );
	return newRoot->template Equilibrate< TSelf >( tmp );
}

template< class TValue >
template< class TSelf >
TSelf *AVLTreeNode< TValue >
::Remove( TSelf *X )
{
	if( X == NULL )
	{
		return dynamic_cast< TSelf* >( this );
	}
	if( this == NULL )
	{
		return NULL;
	}
	TSelf *Y, *Z;
	if( X->GetChild( Left ) != NULL && X->GetChild( Right ) != NULL )
	{
		if( X->GetChild( Left ) != NULL )
		{
			Y = X->template GetChild< TSelf >( Left )->template GetUppestNode< TSelf >();
		}
		else
		{
			Y = X->template GetChild< TSelf >( Right )->template GetLowestNode< TSelf >();
		}
		X->Set( *Y );
	}
	else
	{
		Y = X;
	}
	TSelf *parent = Y->template GetParent< TSelf >();
	TSelf* newRoot = Superclass::template Remove< TSelf >( Y );
	if( parent == NULL )
	{
		parent = this;
	}
	return newRoot->template Equilibrate< TSelf >( parent, true );
}

template< class TValue >
template< class TSelf >
TSelf* AVLTreeNode< TValue >
::Merge( TSelf *root )
{
	if( root == NULL )
	{
		return dynamic_cast< TSelf* >( this );
	}
	if( root->GetParent() != NULL )
	{
		root->GetParent()->AddChild( root->GetParent()->ChildPosition( root ), NULL );
		root->SetParent( NULL );
	}
	if( this == NULL )
	{
		std::cout << " ici " << std::endl;
		return root;
	}
	TSelf *newRoot = dynamic_cast< TSelf* >( this );
	if( root->GetDepth() > this->GetDepth() )
	{
		TSelf *tmp = dynamic_cast< TSelf* >( this );
		newRoot = root;
		newRoot->SetParent( NULL );
		root = tmp;
	}
	TSelf *n, *r;
	bool reverse = *root < *newRoot;
	if( reverse )
	{
		n = root->template GetUppestNode< TSelf >();
	}
	else
	{
		n = root->template GetLowestNode< TSelf >();
	}
	Self *np = n->template GetParent< TSelf >();
	if( reverse )
	{
		r = newRoot->template GetLowestNode< TSelf >( root->GetDepth() );
	}
	else
	{
		r = newRoot->template GetUppestNode< TSelf >( root->GetDepth() );
	}
	if( np != NULL )
	{
		np->AddChild( np->ChildPosition( n ), NULL );
	}
	else
	{
		root = root->template GetChild< TSelf >( Left );
	}
	if( r != NULL  && r->GetParent() != NULL )
	{
		r->GetParent()->AddChild( r->GetParent()->ChildPosition( r ), n );
	}
	if( reverse )
	{
		n->AddChild( Left, root );
		n->AddChild( Right, r );
	}
	else
	{
		n->AddChild( Right, root );
		n->AddChild( Left, r );
	}
	std::cout << " ici " << std::endl;
	if( !n->HasChildren() )
	{
		n->UpdateDepthToRoot( newRoot );
		return newRoot->Equilibrate( n, true );
	}
	std::cout << " ici " << std::endl;
	if( n->GetChild( Left ) != NULL )
	{
		n->GetChild( Left )->UpdateDepthToRoot( newRoot );
	}
	else
	{
		n->GetChild( Right )->UpdateDepthToRoot( newRoot );
	}
	std::cout << " ici " << std::endl;
	return newRoot->template Equilibrate< TSelf >( n, true );
}

template< class TValue >
template< class TSelf >
TSelf* AVLTreeNode< TValue >
::Equilibrate( TSelf *node, bool recursive )
{
	std::cout << " eq " << node << std::endl;
	if( node == NULL )
	{
		std::cout << " ici " << std::endl;
		return dynamic_cast< TSelf* >( this );
	}
	if( this == NULL )
	{
		return dynamic_cast< TSelf* >( this );
	}
	TSelf *L = node->template GetChild< TSelf >( Left ), *R = node->template GetChild< TSelf >( Right );
	TSelf *newRoot = dynamic_cast< TSelf* >( this );
	int balance = L->GetDepth() - R->GetDepth();
	std::cout << " balance " << balance << std::endl;
	std::cout << " L " << L->GetDepth() << " R " << R->GetDepth() << std::endl;
	std::cout << " node " << node->GetDepth() << std::endl;
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

}// namespace itk

#endif
