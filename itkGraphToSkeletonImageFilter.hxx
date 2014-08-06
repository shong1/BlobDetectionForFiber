#ifndef __itkGraphToSkeletonImageFilter_hxx
#define __itkGraphToSkeletonImageFilter_hxx

#include "itkGraphToSkeletonImageFilter.h"

namespace itk
{

template< class TInputGraph, class TOutputImage >
GraphToSkeletonImageFilter< TInputGraph, TOutputImage >
::GraphToSkeletonImageFilter()
{
	this->ProcessObject::SetNumberOfRequiredOutputs( 1 );

	OutputImagePointer output0 = dynamic_cast< OutputImageType * > ( this->MakeOutput( 0 ).GetPointer() );
	this->ProcessObject::SetNthOutput( 0, output0.GetPointer() );
}

template< class TInputGraph, class TOutputImage >
DataObject::Pointer GraphToSkeletonImageFilter< TInputGraph, TOutputImage >
::MakeOutput( unsigned int idx )
{
	OutputImagePointer output = OutputImageType::New();
	return dynamic_cast< DataObject * > ( output.GetPointer() );
}

template< class TInputGraph, class TOutputImage >
void GraphToSkeletonImageFilter< TInputGraph, TOutputImage >
::GenerateOutputInformation()
{
	Superclass::GenerateOutputInformation();
	OutputImagePointer output = OutputImageType::New();
	output->SetRegions( m_Region );
	output->Allocate();
	output->FillBuffer( NumericTraits< PixelType >::Zero );
	output->Update();
	this->GetOutput(0)->Graft( output );

}

template< class TInputGraph, class TOutputImage >
void GraphToSkeletonImageFilter< TInputGraph, TOutputImage >
::GenerateData()
{

	for( typename NodeContainerType::ConstIterator it = this->GetInput()->GetNodeContainer()->Begin();
			it != this->GetInput()->GetNodeContainer()->End();
			++it )
	{
		NodeType node = it.Value();
		for( typename GraphTraitsType::IndexContainerType::iterator idxit = node.LineIndexes.begin();
				idxit != node.LineIndexes.end();
				++idxit )
		{
			if( it.Value().LineIndexes.size() )
			{
				this->GetOutput()->SetPixel( *idxit, it.Value().Weight );
			}
		}
	}
}

} // namespace itk

#endif
