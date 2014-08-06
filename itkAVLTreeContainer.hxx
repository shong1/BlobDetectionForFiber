#ifndef __itkAVLTreeContainer_hxx
#define __itkAVLTreeContainer_hxx

#include "itkAVLTreeContainer.h"

namespace itk
{

template< class TTreeNode >
typename AVLTreeContainer< TTreeNode >
::TreeNodeType * AVLTreeContainer< TTreeNode >
::Add( const ValueType element )
{
	TreeNodeType* node = dynamic_cast< TreeNodeType * >( Superclass::Add( element ) );
	this->Equilibrate( node );
	return node;
}

template< class TTreeNode >
void AVLTreeContainer< TTreeNode >
::Add( TreeNodeType * node )
{
	Superclass::Add( node );
	this->Equilibrate( node );
}

template< class TTreeNode >
void AVLTreeContainer< TTreeNode >
::Merge( Self* tree )
{
	this->Merge( tree->GetRoot() );
	tree->SetRoot( nullptr );
	/*if( tree->GetRoot()->Get() > this->GetRoot()->Get() )
	{
		TreeNodeBaseType *tmp = this->GetRoot();
		this->SetRoot( tree->GetRoot() );
		tree->SetRoot( tmp );
	}
	// http//:stackoverflow.com/questions/2037212/concatenating-merging-joining-two-avl-trees
	TreeNodeType *n, *r;
	n = tree->GetUppestNode();
	tree->Remove( n );
	r = this->GetLowestNode( tree->Depth() );
	r->GetParent()->AddChild( r->GetParent()->ChildPosition( r ), n );
	n->AddChild( Left, dynamic_cast< TreeNodeType* >( tree->GetRoot() ) );
	n->AddChild( Right, r );
	tree->SetRoot( nullptr );
	n->GetChild( Left )->UpdateDepthToRoot();
	n->GetChild( Right )->UpdateDepthToRoot();
	this->Equilibrate( n, true );
*/}

template< class TTreeNode >
void AVLTreeContainer< TTreeNode >
::Merge( TreeNodeType * root )
{
	if( root == NULL )
	{
		return;
	}
	if( root->GetParent() != NULL )
	{
		root->SetParent( NULL );
	}
	if( this->GetRoot() == NULL )
	{
		this->SetRoot( root );
		return;
	}
	if( root->GetDepth() > this->GetRoot()->GetDepth() )
	{
		TreeNodeType * tmp = this->GetRoot();
		this->SetRoot( root );
		root = tmp;
		if( root == NULL )
		{
			return;
		}
	}
	TreeNodeType *n, *r;
	bool reverse = root->Get() < this->GetRoot()->Get();
	if( reverse )
	{
		n = root->GetUppestNode();
	}
	else
	{
		n = root->GetLowestNode();
	}
	TreeNodeType *np = n->GetParent();
	if( reverse )
	{
		r = this->GetLowestNode( root->GetDepth() );
	}
	else
	{
		r = this->GetUppestNode( root->GetDepth() );
	}
	if( np != NULL )
	{
		np->AddChild( np->ChildPosition( n ), NULL );
	}
	else
	{
		root = root->GetChild( Left );
	}
	if( r != NULL && r->GetParent() != NULL )
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
	if( n->GetChild( Left ) == NULL && n->GetChild( Right ) == NULL )
	{
		n->UpdateDepthToRoot( this->GetRoot() );
		this->Equilibrate( n, true );
		return;
	}
	if( n->GetChild( Left ) != NULL )
	{
		n->GetChild( Left )->UpdateDepthToRoot( this->GetRoot() );
	}
	if( n->GetChild( Right ) != NULL )
	{
		n->GetChild( Right )->UpdateDepthToRoot( this->GetRoot() );
	}
	this->Equilibrate( n, true );
}

template< class TTreeNode >
typename AVLTreeContainer< TTreeNode >
::Self* AVLTreeContainer< TTreeNode >
::SizeSplitLeft( int size, Self *tree )
{
	if( this->GetRoot() == NULL )
	{
		return NULL;
	}
	TreeNodeType *node = this->GetRoot();
	this->SetRoot( nullptr );
	TreeNodeType *tmp;
	int lsize = 0;
	while( node != NULL )
	{
		if( node->GetChild( Left )->GetSize() + lsize < size )
		{
			tmp = node->GetChild( Right );
			lsize += node->GetChild( Left )->GetSize() + 1;
			tree->Merge( node->GetChild( Left ) );
			tree->Add( node );
		}
		else
		{
			tmp = node->GetChild( Left );
			this->Merge( node->GetChild( Right ) );
			this->Add( node );
		}
		node = tmp;
	}
	return tree;
}

template< class TTreeNode >
typename AVLTreeContainer< TTreeNode >
::Self* AVLTreeContainer< TTreeNode >
::SizeSplitRight( int size, Self *tree )
{
	if( this->GetRoot() == NULL )
	{
		return NULL;
	}
	TreeNodeType *node = this->GetRoot();
	this->SetRoot( nullptr );
	TreeNodeType *tmp;
	int lsize = 0;
	while( node != NULL )
	{
		if( node->GetChild( Right )->GetSize() + lsize < size )
		{
			tmp = node->GetChild( Left );
			lsize += node->GetChild( Right )->GetSize() + 1;
			tree->Merge( node->GetChild( Right ) );
			tree->Add( node );
		}
		else
		{
			tmp = node->GetChild( Right );
			this->Merge( node->GetChild( Left ) );
			this->Add( node );
		}
		node = tmp;
	}
	return tree;
}

template< class TTreeNode >
typename AVLTreeContainer< TTreeNode >
::Self* AVLTreeContainer< TTreeNode >
::Split( const ValueType element, Self * tree )
{
	std::cout << " Split " << std::endl;
	// this right, tree left
	if( this->GetRoot() == NULL )
	{
		return NULL;
	}
	TreeNodeType *node = this->GetRoot();
	this->SetRoot( nullptr );
//	tree->SetRoot( nullptr );
	
	TreeNodeType *tmp;
	while( node != NULL )
	{
		if( node->Get() < element )
		{
			tmp = node->GetChild( Right );
			tree->Merge( node->GetChild( Left ) );
			tree->Add( node );
		}
		else
		{
			tmp = node->GetChild( Left );
			this->Merge( node->GetChild( Right ) );
			this->Add( node );
		}
		node = tmp;
	}

	return tree;
/*

	Pointer output = Self::New();
	TreeNodeType * tmp = this->GetRoot();
	ChildIdentifier side = Left;
	if( node->Get() > this->GetRoot()->Get() )
	{
		side = Right;
		while( node->Get() > tmp->Get() )
		{
			tmp = tmp->GetChild( Right );
		}
	}
	else
	{
		while( node->Get() < tmp->Get() )
		{
			tmp = tmp->GetChild( Left );
		}
	}
	TreeNodeType * subtreeRoot = tmp->GetChild( !side );
	tmp->AddChild( !side, NULL );
	output->SetRoot( tmp->GetChild( side ) );
	output->Add( tmp );

	int balance = std::abs( this->GetRoot()->GetDepth() - subtreeRoot->GetDepth() );

	while( balance > 1 )
	{
		this->Rotate( this->GetRoot(), side );
	}

	return output;
*/}

template< class TTreeNode >
bool AVLTreeContainer< TTreeNode >
::Remove( const ValueType element )
{
	return this->Remove( this->GetNode( element ) );
}

template< class TTreeNode >
bool AVLTreeContainer< TTreeNode >
::Remove( TreeNodeType * X )
{
	if( X == NULL )
	{
		return false;
	}
	TreeNodeType * Y, * Z;
	if( X->GetChild( Left ) != NULL && X->GetChild( Right ) != NULL )
	{
		if( X->GetChild( Left ) != NULL )
		{
			Y = this->GetUppestNode( X->GetChild( Left ) );
		}
		else
		{
			Y = this->GetLowestNode( X->GetChild( Right ) );
		}
		X->Set( Y->Get() );
	}
	else
	{
		Y = X;
	}
	TreeNodeType * parent = Y->GetParent();
	if( !Superclass::Remove( Y ) )
	{
		return false;
	}
	if( parent == NULL )
	{
		parent = Superclass::m_Root;
	}
	return this->Equilibrate( parent, true );
}

template< class TTreeNode >
bool AVLTreeContainer< TTreeNode >
::Equilibrate( TreeNodeType * node, bool recursive )
{
	return this->Equilibrate( node, this->GetRoot(), recursive );
}

template< class TTreeNode >
bool AVLTreeContainer< TTreeNode >
::Equilibrate( TreeNodeType * node, TreeNodeType * root, bool recursive )
{
	if( node == NULL || root != NULL && node == root->GetParent() )
	{
		return true;
	}
	TreeNodeType *L = node->GetChild( Left ), *R = node->GetChild( Right );
	int balance = this->Depth( L ) - this->Depth( R );
	if( balance == 2 )
	{
		TreeNodeType *LL = L->GetChild( Left ), *LR = L->GetChild( Right );
		if( this->Depth( LL ) > this->Depth( LR ) )
		{
			if( !this->RotateRight( node ) )
			{
				return false;
			}
		}
		else
		{
			if( !this->RotateLeft( L ) )
			{
				return false;
			}
			if( !this->RotateRight( node ) )
			{
				return false;
			}
		}
		if( !recursive )
		{
			return true;
		}
	}
	if( balance == -2 )
	{
		TreeNodeType *RL = R->GetChild( Left ), *RR = R->GetChild( Right );
		if( this->Depth( RR ) > this->Depth( RL ) )
		{
			if( !this->RotateLeft( node ) )
			{
				return false;
			}
		}
		else
		{
			if( !this->RotateRight( R ) )
			{
				return false;
			}
			if( !this->RotateLeft( node ) )
			{
				return false;
			}
		
		}
		if( !recursive )
		{
			return true;
		}
	}
	this->Equilibrate( node->GetParent(), root, recursive );
}

}// namespace itk

#endif
