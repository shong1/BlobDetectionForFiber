// THIS PROJECT USES C++11
// NO LOWER VERSION
//
#include <iostream>
#include <fstream>
#include <math.h>
#include <list>
#include <stack>
#include <string>


using namespace std;
using std::ofstream;
using std::ifstream;


	/*
	 * Image
	 */

#include "itkImage.h"

const unsigned int Dim = 3;
typedef float RealType;
typedef unsigned char BinType;



typedef unsigned short 					PixelType;
typedef itk::Image< PixelType, Dim > 	ImageType;
typedef itk::Image< RealType, Dim > 	RealImageType;
typedef itk::ImageRegion< Dim >			ImageRegionType;
typedef itk::Image< BinType, Dim > 		LabelType;

typedef itk::Image< unsigned long , Dim > 	SegLabelType;


const PixelType MaxPixelValue = itk::NumericTraits< PixelType >::max();
const PixelType MinPixelValue = itk::NumericTraits< PixelType >::min();

#include "itkSkeletonGraph.h"
#include "itkLineGraphTraits.h"

typedef itk::LineGraphTraits< Dim >				GraphTraitsType;
typedef itk::SkeletonGraph< GraphTraitsType >	GraphType;

	/*
	 * Inputs and Outputs
	 */

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

typedef itk::ImageFileReader< ImageType > ReaderType;
typedef itk::ImageFileReader< RealImageType > RealReaderType;
typedef itk::ImageFileWriter< ImageType > WriterType;
typedef itk::ImageFileWriter< RealImageType >							RealWriterType;
typedef itk::ImageFileWriter< LabelType > LabelWriterType;

typedef itk::ImageFileWriter<  SegLabelType > SegLabelWriterType;

	/*
	 * Filters
	 */

#include "itkBinaryThinningImageFilter3D.h"
#include "itkSkeletonToGraphFilter.h"
#include "itkGraphToSkeletonImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkFiberExtractionGraphFilter.h"
#include "itkLineShapeImageFilter.h"

#include "itkAVLIntervalTreeNode.h"
#include "itkAVLVectorIntervalTreeNode.h"
#include "itkHalfConvexHullTreeNode.h"
#include "itkHalfConvexHullTreeContainer.h"
#include "itkConvexHull.h"
#include "itkBlurred3DLineSegmentation.hxx"
#include "itkImageDuplicator.h"

#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkScalarImageKmeansImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkStatisticsImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkGrayscaleErodeImageFilter.h"
#include "itkGrayscaleDilateImageFilter.h"


#include "itkGradientMagnitudeImageFilter.h"
#include "itkWatershedImageFilter.h"

typedef itk::Vector< int, 2 >											VectorType;
typedef itk::AVLIntervalTreeNode< VectorType >									AVLTreeNodeType;
typedef itk::HalfConvexHullTreeNode< VectorType >					HalfConvexHullTreeNodeType;
typedef itk::HalfConvexHullTreeContainer< HalfConvexHullTreeNodeType >			HalfConvexHullTreeContainerType;
typedef itk::BinaryTreeContainer< AVLTreeNodeType >						BinaryTreeContainerType;
typedef AVLTreeNodeType													TreeNodeType;

typedef itk::LineShapeImageFilter< ImageType > 							LineShapeFilterType;
typedef itk::BinaryThinningImageFilter3D< ImageType, ImageType >		BinaryThinningFilterType;


typedef itk::SkeletonToGraphFilter< ImageType, GraphType, RealImageType, ImageType >	SkeletonToGraphFilterType;
typedef itk::GraphToSkeletonImageFilter< GraphType, RealImageType >			GraphToSkeletonFilterType;
typedef itk::FiberExtractionGraphFilter< GraphType >					FiberExtractionFilterType;

typedef itk::MultiplyImageFilter< ImageType, RealImageType, RealImageType > MultiplyFilterType;

typedef itk::ImageDuplicator< RealImageType >	DuplicatorType;

typedef itk::MinimumMaximumImageCalculator< ImageType >					MinMaxFilterType;
typedef itk::StatisticsImageFilter< ImageType > 						StaticFilterType;
typedef itk::ScalarImageKmeansImageFilter< ImageType, LabelType > 		KMeansType;
typedef itk::StatisticsImageFilter< SegLabelType > 						LabelStaticFilterType;
typedef itk::StatisticsImageFilter< RealImageType > 							RealStatisticFilterType;

typedef itk::ConnectedComponentImageFilter< LabelType, ImageType >		CCAFilterType;
typedef itk::RelabelComponentImageFilter< ImageType, LabelType > 		RelabelFilterType;
typedef itk::RelabelComponentImageFilter< ImageType, SegLabelType > 		SegRelabelFilterType;

typedef itk::BinaryThresholdImageFilter< ImageType, LabelType > 		BinaryFilterType;
typedef itk::BinaryBallStructuringElement< LabelType::PixelType, 3 >	StructureElementType;
typedef itk::BinaryErodeImageFilter< LabelType, LabelType, StructureElementType > 			BinErodeFilterType;
typedef itk::GrayscaleDilateImageFilter< ImageType, ImageType, StructureElementType >		GrayDilateFilterType;
typedef itk::SubtractImageFilter< ImageType, ImageType > 				SubtractorType;

typedef itk::GradientMagnitudeImageFilter< ImageType, RealImageType > 	GradientMagnitudeFilterType;
typedef itk::WatershedImageFilter< RealImageType > 						WatershedFilterType;

#include "itkN4BiasFieldCorrectionImageFilter.h"
typedef itk::N4BiasFieldCorrectionImageFilter< ImageType, LabelType, ImageType > 	N4CorrectorType;

#include "vtkSmartPointer.h"
#include "vtkDoubleArray.h"
#include "vtkKMeansStatistics.h"
#include "vtkIntArray.h"
#include "vtkTable.h"
#include "vtkVariant.h"

//#include "vtkImageData.h"
//#include "vtkPointData.h"
//#include "vtkFloatArray.h"
//#include "vtkPolyData.h"
//#include "vtkUnsignedShortArray.h"
//#include "vtkDataObject.h"
//#include "vtkContourFilter.h"
//#include "vtkLookupTable.h"
//#include "vtkImageMapToColors.h"
//#include "vtkProbeFilter.h"
//#include "vtkSmoothPolyDataFilter.h"
//#include "vtkPolyDataMapper.h"
//#include "vtkScalarBarActor.h"
//#include "vtkTextProperty.h"
//#include "vtkTextMapper.h"
//#include "vtkActor2D.h"
//#include "vtkActor.h"
//#include "vtkRenderer.h"
//#include "vtkRenderWindow.h"
//#include "vtkRenderWindowInteractor.h"
//#include "vtkStructuredPointsReader.h"
//

#define VTK_CREATE( type, var ) vtkSmartPointer< type > var = vtkSmartPointer< type >::New()

void histoGraphMap( GraphType::Pointer graph );

void display( TreeNodeType* node );

int main( int argc, char* argv[] )
{
	/*
	 * Checking input parametersmain
	 * The parameter is the file without the extension
	 */

//	if( argc < 2 )
//	{
//		std::cout << "Error : Not enough input arguments" << std::endl;
//		return EXIT_FAILURE;
//	}
//
//	std::string str = std::string( argv[1] );

	std::string str = "/home/shong/Research/FibersInConcrete/Yohan/data/FlexFibers_cropped";
	std::string strRes = "/home/shong/Research/FibersInConcrete/Yohan/data/FlexFibers_cropped_Res/FlexFibers_cropped";

	WriterType::Pointer writer = WriterType::New();
	RealWriterType::Pointer realWriter = RealWriterType::New();
	LabelWriterType::Pointer labelWriter = LabelWriterType::New();
	SegLabelWriterType::Pointer segLabelWriter = SegLabelWriterType::New();

	/*
	 *	Reading input image
	 */
/*
	// The file format is nrrd, but it can be changed
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( str + ".nrrd" );
	reader->Update();

	ImageType::Pointer input = reader->GetOutput();

	StaticFilterType::Pointer statFilter = StaticFilterType::New();
	statFilter->SetInput( input );
	statFilter->Update();
	float meanI = statFilter->GetMean();
	float stdI = statFilter->GetSigma();
	float minImg = statFilter->GetMinimum();
	float maxImg = statFilter->GetMaximum();
	cout << "Mean+ : " << meanI << " Std : " << stdI << " Min : " << minImg << " Max : " << maxImg << endl;

	float statMin = meanI - 6 * stdI;
	float statMax = meanI + 6 * stdI;
	cout << "StatMin : " << statMin << " StatMax : " << statMax << endl;

	ImageType::RegionType region = input->GetLargestPossibleRegion();
	ImageType::SizeType dim3 = region.GetSize();

	for( int x = 0; x < dim3[0]; x++ )
	for( int y = 0; y < dim3[1]; y++ )
	for( int z = 0; z < dim3[2]; z++ )
	{
		ImageType::IndexType idx;
		idx[0] = x;
		idx[1] = y;
		idx[2] = z;

		unsigned short val2 = input->GetPixel( idx );

		if( val2 < statMin )
			input->SetPixel( idx, statMin );
		else if( val2 > statMax )
			input->SetPixel( idx, statMax );
	}

	input->Update();

//	writer->SetInput( input );
//	writer->SetFileName( str + "_IntensityCut.nrrd" );
//	writer->Update();
*/


/////////////////////////////////////////////////////////////////////
// 					Image Correction 							   //
/////////////////////////////////////////////////////////////////////
/*
	std::cout << " ***** Start " << std::endl;

	// N3 MRI Bias Field
	LabelType::Pointer N4Label = LabelType::New();
	N4Label->SetLargestPossibleRegion( region );
	N4Label->Allocate();

	for( int x = 0; x < dim3[0]; x++ )
	for( int y = 0; y < dim3[1]; y++ )
	for( int z = 0; z < dim3[2]; z++ )
	{
		LabelType::IndexType idx;
		idx[0] = x;
		idx[1] = y;
		idx[2] = z;

		N4Label->SetPixel( idx, 1 );
	}
	N4Label->Update();

	std::cout << " ***** Image Correction " << std::endl;

	N4CorrectorType::Pointer corrector = N4CorrectorType::New();
	corrector->SetInput( input );
	corrector->SetMaskImage( N4Label );
	corrector->Update();

	ImageType::Pointer N4CorrectedImg = corrector->GetOutput();
	writer->SetInput( N4CorrectedImg );
	writer->SetFileName( str + "_N4Corrected.nrrd" );
	writer->Update();

	input = corrector->GetOutput();
*/
/*
	ReaderType::Pointer readerCorrected = ReaderType::New();
	readerCorrected->SetFileName( str + "_N4Corrected.nrrd" );
	readerCorrected->Update();
	ImageType::Pointer input = readerCorrected->GetOutput();

	ImageType::RegionType region = input->GetLargestPossibleRegion();
	ImageType::SizeType dim3 = region.GetSize();

	StaticFilterType::Pointer statFilter = StaticFilterType::New();
	statFilter->SetInput( input );
	statFilter->Update();
	float meanI = statFilter->GetMean();
	float stdI = statFilter->GetSigma();
	float minImg = statFilter->GetMinimum();
	float maxImg = statFilter->GetMaximum();
	cout << "Mean+ : " << meanI << " Std : " << stdI << " Min : " << minImg << " Max : " << maxImg << endl;

	float statMin = meanI - 6 * stdI;
	float statMax = meanI + 6 * stdI;


/////////////////////////////////////////////////////////////////////
//				Extract Image Information 						  ///
/////////////////////////////////////////////////////////////////////
	// Watershed Segmentation for Blob Detection
	GradientMagnitudeFilterType::Pointer gradientFilter = GradientMagnitudeFilterType::New();
	gradientFilter->SetInput( input );
	gradientFilter->Update();

	realWriter->SetInput( gradientFilter->GetOutput() );
	realWriter->SetFileName( str + "_gradMag.nrrd" );
	realWriter->Update();
	RealImageType::SizeType dim = region.GetSize();

	RealImageType::Pointer vesselImgMultiSigma = RealImageType::New();
	vesselImgMultiSigma->SetRegions( input->GetLargestPossibleRegion() );
	vesselImgMultiSigma->SetSpacing( input->GetSpacing() );
	vesselImgMultiSigma->SetOrigin( input->GetOrigin() );
	vesselImgMultiSigma->Allocate();

	RealImageType::Pointer vesselImgMultiSigmaSato = RealImageType::New();
	vesselImgMultiSigmaSato->SetRegions( input->GetLargestPossibleRegion() );
	vesselImgMultiSigmaSato->SetSpacing( input->GetSpacing() );
	vesselImgMultiSigmaSato->SetOrigin( input->GetOrigin() );
	vesselImgMultiSigmaSato->Allocate();

	RealImageType::Pointer blobImgMultiSigma = RealImageType::New();
	blobImgMultiSigma->SetRegions( input->GetLargestPossibleRegion() );
	blobImgMultiSigma->SetSpacing( input->GetSpacing() );
	blobImgMultiSigma->SetOrigin( input->GetOrigin() );
	blobImgMultiSigma->Allocate();

	RealImageType::Pointer eig0ImgMultiSigma = RealImageType::New();
	eig0ImgMultiSigma->SetRegions( input->GetLargestPossibleRegion() );
	eig0ImgMultiSigma->SetSpacing( input->GetSpacing() );
	eig0ImgMultiSigma->SetOrigin( input->GetOrigin() );
	eig0ImgMultiSigma->Allocate();

	RealImageType::Pointer eig1ImgMultiSigma = RealImageType::New();
	eig1ImgMultiSigma->SetRegions( input->GetLargestPossibleRegion() );
	eig1ImgMultiSigma->SetSpacing( input->GetSpacing() );
	eig1ImgMultiSigma->SetOrigin( input->GetOrigin() );
	eig1ImgMultiSigma->Allocate();

	RealImageType::Pointer eig2ImgMultiSigma = RealImageType::New();
	eig2ImgMultiSigma->SetRegions( input->GetLargestPossibleRegion() );
	eig2ImgMultiSigma->SetSpacing( input->GetSpacing() );
	eig2ImgMultiSigma->SetOrigin( input->GetOrigin() );
	eig2ImgMultiSigma->Allocate();


	for( int x = 0; x < dim3[0]; x++ )
	for( int y = 0; y < dim3[1]; y++ )
	for( int z = 0; z < dim3[2]; z++ )
	{
		RealImageType::IndexType idx;
		idx[0] = x;
		idx[1] = y;
		idx[2] = z;

		vesselImgMultiSigma->SetPixel( idx, 0 );
		vesselImgMultiSigmaSato->SetPixel( idx, 0 );
		eig0ImgMultiSigma->SetPixel( idx, 0 );
		eig1ImgMultiSigma->SetPixel( idx, 0 );
		eig2ImgMultiSigma->SetPixel( idx, 0 );
		blobImgMultiSigma->SetPixel( idx, 0 );
	}
	vesselImgMultiSigma->Update();
	vesselImgMultiSigmaSato->Update();
	eig0ImgMultiSigma->Update();
	eig1ImgMultiSigma->Update();
	eig2ImgMultiSigma->Update();
	blobImgMultiSigma->Update();

	std::cout << " ***** Line Shape Filter " << std::endl;
	// This filter does a hessian + gaussian filter
	// We have then access to eigen values
	// Scale Space Search
	for( int i = 5; i < 26; i++ )
	{
		double sigma = i * 0.2;
		string strS = str + std::to_string( sigma );

		LineShapeFilterType::Pointer lineShapeFilter = LineShapeFilterType::New();
		lineShapeFilter->SetInput( input );
		// this parameter has to be about the width of a fiber
		lineShapeFilter->SetSigma( sigma );
		lineShapeFilter->SetExtractBrightLine( false );
		lineShapeFilter->EigenValuesExtractionOn();
		lineShapeFilter->LabelImageOn();
		lineShapeFilter->Update();

//		// Binary image of the fiber-like shapesppt
//		writer->SetInput( lineShapeFilter->GetBinaryOutput() );
//		writer->SetFileName( strS + "_binary.nrrd" );
//		writer->Update();
//
//		// Label image depending of the sign of the 3 sorted eigenvalues
//		writer->SetInput( lineShapeFilter->GetLabelOutput() );
//		writer->SetFileName( strS + "_label.nrrd" );
//		writer->Update();

		RealImageType::Pointer eigImg0 = lineShapeFilter->GetEigenValuesOutput( 0 );
		RealImageType::Pointer eigImg1 = lineShapeFilter->GetEigenValuesOutput( 1 );
		RealImageType::Pointer eigImg2 = lineShapeFilter->GetEigenValuesOutput( 2 );

		RealStatisticFilterType::Pointer eigStatFilter = RealStatisticFilterType::New();
		eigStatFilter->SetInput( eigImg0 );
		eigStatFilter->Update();
		float minEig0 = eigStatFilter->GetMinimum();
		float maxEig0 = eigStatFilter->GetMaximum();

		cout << " Min : " << minEig0 << " Max : " << maxEig0 << endl;

		eigStatFilter->SetInput( eigImg1);
		eigStatFilter->Update();
		float minEig1 = eigStatFilter->GetMinimum();
		float maxEig1 = eigStatFilter->GetMaximum();

		cout << " Min : " << minEig1 << " Max : " << maxEig1 << endl;

		eigStatFilter->SetInput( eigImg2 );
		eigStatFilter->Update();
		float minEig2 = eigStatFilter->GetMinimum();
		float maxEig2 = eigStatFilter->GetMaximum();

		cout << " Min : " << minEig2 << " Max : " << maxEig2 << endl;

		double rEig0Inten = ( statMax - statMin ) / ( maxEig0 - minEig0 );
		double rEig1Inten = ( statMax - statMin ) / ( maxEig1 - minEig1 );
		double rEig2Inten = ( statMax - statMin ) / ( maxEig2 - minEig2 );

//		realWriter->SetFileName( strS + "_eig0.nrrd" );
//		realWriter->SetInput( eigImg0 );
//		realWriter->Update();
//
//		realWriter->SetFileName( strS + "_eig1.nrrd" );
//		realWriter->SetInput( eigImg1 );
//		realWriter->Update();
//
//		realWriter->SetFileName( strS + "_eig2.nrrd" );
//		realWriter->SetInput( eigImg2 );
//		realWriter->Update();
//
//		// Vesselness output (fiber-shape likelihood)
//		realWriter->SetInput( lineShapeFilter->GetVesselnessOutput() );
//		realWriter->SetFileName( strS + "_vesselness.nrrd" );
//		realWriter->Update();
//
//		// Label image depending of the sign of the 3 sorted eigenvalues
//		writer->SetInput( lineShapeFilter->GetLabelOutput() );
//		writer->SetFileName( strS + "_label.nrrd" );
//		writer->Update();

		RealImageType::Pointer blobImg = RealImageType::New();
		blobImg->SetRegions( eigImg0->GetLargestPossibleRegion() );
		blobImg->SetSpacing( eigImg0->GetSpacing() );
		blobImg->SetOrigin( eigImg0->GetOrigin() );
		blobImg->Allocate();

		LabelType::Pointer blobBin = LabelType::New();
		blobBin->SetRegions( eigImg0->GetLargestPossibleRegion() );
		blobBin->SetSpacing( eigImg0->GetSpacing() );
		blobBin->SetOrigin( eigImg0->GetOrigin() );
		blobBin->Allocate();

		LabelType::Pointer vesselBin = LabelType::New();
		vesselBin->SetRegions( eigImg0->GetLargestPossibleRegion() );
		vesselBin->SetSpacing( eigImg0->GetSpacing() );
		vesselBin->SetOrigin( eigImg0->GetOrigin() );
		vesselBin->Allocate();

		RealImageType::Pointer vesselImg = RealImageType::New();
		vesselImg->SetRegions( eigImg0->GetLargestPossibleRegion() );
		vesselImg->SetSpacing( eigImg0->GetSpacing() );
		vesselImg->SetOrigin( eigImg0->GetOrigin() );
		vesselImg->Allocate();

		cout << dim3[0] << ", " << dim3[1] << ", " << dim3[2] << endl;
		float alpha = 0.2;
		float beta = 0.2;
		float c = sqrt( maxEig0 * maxEig0 + maxEig1 * maxEig1 + maxEig2 * maxEig2 ) / 2.0;

		for( int x = 0; x < dim3[0]; x++ )
		for( int y = 0; y < dim3[1]; y++ )
		for( int z = 0; z < dim3[2]; z++ )
		{
			RealImageType::IndexType idx;
			idx[0] = x;
			idx[1] = y;
			idx[2] = z;

			float lambda1 = eigImg0->GetPixel( idx ) * sigma * sigma;
			float lambda2 = eigImg1->GetPixel( idx ) * sigma * sigma;
			float lambda3 = eigImg2->GetPixel( idx ) * sigma * sigma;

			// eigImg2 > eigImg1 - 1.0 > eigImg0
			// lambda3 > lambda2 > lambda1
			float blVal = abs( lambda1 ) / sqrt( abs( lambda2 * lambda3 ) );
			float sVal = sqrt( lambda1 * lambda1 + lambda2 * lambda2 + lambda3 * lambda3 );
			float aVal = abs( lambda2 / lambda3 );

			float V0 = ( 1 - exp( -( blVal * blVal ) / ( 2 * beta * beta ) ) ) * ( 1 - exp( -( sVal * sVal ) / ( 2 * c * c ) ) );

			float Ves = ( 1 - exp( -( ( aVal ) * ( aVal ) ) / ( 2 * alpha * alpha ) ) ) * ( exp( -( blVal * blVal ) / ( 2 * beta * beta ) ) ) * ( 1 - exp( -( sVal * sVal ) / ( 2 * c * c ) ) );
			Ves = Ves * 10000;
			V0 = V0 * 10000;

			//		float V0 = 1 - exp( ( -blVal * blVal ) / ( 2 * beta * beta ) );
			blobImg->SetPixel( idx, V0 );
			vesselImg->SetPixel( idx, Ves );

			double VesSato = lineShapeFilter->GetVesselnessOutput()->GetPixel( idx );

			if( V0 > blobImgMultiSigma->GetPixel( idx ) )
				blobImgMultiSigma->SetPixel( idx, V0 );

			if( Ves > vesselImgMultiSigma->GetPixel( idx ) )
				vesselImgMultiSigma->SetPixel( idx, Ves );

			if( VesSato > vesselImgMultiSigmaSato->GetPixel( idx ) )
				vesselImgMultiSigmaSato->SetPixel( idx, VesSato );

			if( abs( lambda1 ) > abs( eig0ImgMultiSigma->GetPixel( idx ) ) )
				eig0ImgMultiSigma->SetPixel( idx, lambda1 );

			if( abs( lambda2 ) > abs( eig1ImgMultiSigma->GetPixel( idx ) ) )
				eig1ImgMultiSigma->SetPixel( idx, lambda2 );

			if( abs( lambda3 ) > abs( eig2ImgMultiSigma->GetPixel( idx ) ) )
				eig2ImgMultiSigma->SetPixel( idx, lambda3 );

			if( V0 > 0.2 && lambda1 < 0 && lambda2 < 0 && lambda3 < 0 )
				blobBin->SetPixel( idx, 255 );
			else
				blobBin->SetPixel( idx, 0 );

			int vVal = lineShapeFilter->GetLabelOutput()->GetPixel( idx );

			if( vVal == 0 )
				vesselBin->SetPixel( idx, 1 );
			else
				vesselBin->SetPixel( idx, 0 );
		}

		vesselBin->Update();
		blobImg->Update();
		blobBin->Update();
		vesselImg->Update();

//		realWriter->SetInput( blobImg );
//		realWriter->SetFileName( strS + "blobness.nrrd" );
//		realWriter->Update();
//
//		realWriter->SetInput( vesselImg );
//		realWriter->SetFileName( strS + "vesselnessF.nrrd" );
//		realWriter->Update();
//
//		labelWriter->SetInput( vesselBin );
//		labelWriter->SetFileName( strS+"_labelV.nrrd" );
//		labelWriter->Update();
	}

	vesselImgMultiSigmaSato->Update();
	vesselImgMultiSigma->Update();
	eig0ImgMultiSigma->Update();
	eig1ImgMultiSigma->Update();
	eig2ImgMultiSigma->Update();
	blobImgMultiSigma->Update();


	realWriter->SetInput( vesselImgMultiSigma );
	realWriter->SetFileName( str + "_vesselnessF_Scale.nrrd" );
	realWriter->Update();

	realWriter->SetInput( vesselImgMultiSigmaSato );
	realWriter->SetFileName( str + "_vesselnessSato_Scale.nrrd" );
	realWriter->Update();

	realWriter->SetInput( blobImgMultiSigma );
	realWriter->SetFileName( str + "_Blobness_Scale.nrrd" );
	realWriter->Update();

	realWriter->SetInput( eig0ImgMultiSigma );
	realWriter->SetFileName( str + "_eig0_Scale.nrrd" );
	realWriter->Update();
	realWriter->SetInput( eig1ImgMultiSigma );
	realWriter->SetFileName( str + "_eig1_Scale.nrrd" );
	realWriter->Update();
	realWriter->SetInput( eig2ImgMultiSigma );
	realWriter->SetFileName( str + "_eig2_Scale.nrrd" );
	realWriter->Update();
*/

/////////////////////////////////////////////////////////////
// 			Load Preprocessed Image Data 					//
//////////////////////////////////////////////////////////////

// Image Gradient For Watershed Segmentation
	ReaderType::Pointer readerCorrected = ReaderType::New();
	readerCorrected->SetFileName( str + "_N4Corrected.nrrd" );
	readerCorrected->Update();
	ImageType::Pointer input = readerCorrected->GetOutput();

	ImageType::RegionType region = input->GetLargestPossibleRegion();
	ImageType::SizeType dim3 = region.GetSize();

	RealReaderType::Pointer realReaderGrad = RealReaderType::New();
	realReaderGrad->SetFileName( str + "_gradMag.nrrd" );
	realReaderGrad->Update();
	RealImageType::Pointer gradient = realReaderGrad->GetOutput();

	RealReaderType::Pointer realReaderEig0 = RealReaderType::New();
	realReaderEig0->SetFileName( str + "_eig0_Scale.nrrd" );
	realReaderEig0->Update();
	RealImageType::Pointer eigImg0 = realReaderEig0->GetOutput();

	RealReaderType::Pointer realReaderEig1 = RealReaderType::New();
	realReaderEig1->SetFileName( str + "_eig1_Scale.nrrd" );
	realReaderEig1->Update();
	RealImageType::Pointer eigImg1 = realReaderEig1->GetOutput();

	RealReaderType::Pointer realReaderEig2 = RealReaderType::New();
	realReaderEig2->SetFileName( str + "_eig2_Scale.nrrd" );
	realReaderEig2->Update();
	RealImageType::Pointer eigImg2 = realReaderEig2->GetOutput();

	RealReaderType::Pointer realReaderBlob = RealReaderType::New();
	realReaderBlob->SetFileName( str + "_Blobness_Scale.nrrd" );
	realReaderBlob->Update();
	RealImageType::Pointer blobImg = realReaderBlob->GetOutput();

	RealReaderType::Pointer realReaderVessel = RealReaderType::New();
	realReaderVessel->SetFileName( str + "_vesselnessSato_Scale.nrrd" );
	realReaderVessel->Update();
	RealImageType::Pointer vesselImg = realReaderVessel->GetOutput();

	// Image Statistic
	StaticFilterType::Pointer statFilter = StaticFilterType::New();
	statFilter->SetInput( input );
	statFilter->Update();
	float meanI = statFilter->GetMean();
	float stdI = statFilter->GetSigma();
	float minImg = statFilter->GetMinimum();
	float maxImg = statFilter->GetMaximum();
	cout << "Mean+ : " << meanI << " Std : " << stdI << " Min : " << minImg << " Max : " << maxImg << endl;

	float statMin = meanI - 6 * stdI;
	float statMax = meanI + 6 * stdI;

	RealStatisticFilterType::Pointer eigStatFilter = RealStatisticFilterType::New();
	eigStatFilter->SetInput( eigImg0 );
	eigStatFilter->Update();
	float minEig0 = eigStatFilter->GetMinimum();
	float maxEig0 = eigStatFilter->GetMaximum();

	cout << " Min : " << minEig0 << " Max : " << maxEig0 << endl;

	eigStatFilter->SetInput( eigImg1);
	eigStatFilter->Update();
	float minEig1 = eigStatFilter->GetMinimum();
	float maxEig1 = eigStatFilter->GetMaximum();

	cout << " Min : " << minEig1 << " Max : " << maxEig1 << endl;

	eigStatFilter->SetInput( eigImg2 );
	eigStatFilter->Update();
	float minEig2 = eigStatFilter->GetMinimum();
	float maxEig2 = eigStatFilter->GetMaximum();
	cout << " Min : " << minEig2 << " Max : " << maxEig2 << endl;

	double rEig0Inten = ( maxImg - minImg ) / ( maxEig0 - minEig0 );
	double rEig1Inten = ( maxImg - minImg ) / ( maxEig1 - minEig1 );
	double rEig2Inten = ( maxImg - minImg ) / ( maxEig2 - minEig2 );

	eigStatFilter->SetInput( blobImg );
	eigStatFilter->Update();
	float minBlob = eigStatFilter->GetMinimum();
	float maxBlob = eigStatFilter->GetMaximum();
	cout << " Min : " << minBlob << " Max : " << maxBlob << endl;

	eigStatFilter->SetInput( vesselImg );
	eigStatFilter->Update();
	float minVessel = eigStatFilter->GetMinimum();
	float maxVessel = eigStatFilter->GetMaximum();
	cout << " Min : " << minVessel << " Max : " << maxVessel << endl;

	double rVesInten = ( maxImg - minImg ) / ( maxVessel - minVessel );
	double rBlobInten = ( maxImg - minImg ) / ( maxBlob - minBlob );

///////////////////////////////////////////////////////////////////
// Image Processing												 //
///////////////////////////////////////////////////////////////////
	WatershedFilterType::Pointer watersheder = WatershedFilterType::New();
	watersheder->SetInput( gradient );
	watersheder->SetThreshold( 0.001 );
	watersheder->SetLevel( 0.2 );
	watersheder->Update();
	SegLabelType::Pointer waterShedLabel = watersheder->GetOutput();

	segLabelWriter->SetInput( waterShedLabel );
	segLabelWriter->SetFileName( strRes + "_watershed_001_2.nrrd" );
	segLabelWriter->Update();

	// Label Statisitic Filter
	LabelStaticFilterType::Pointer labelStatFilter = LabelStaticFilterType::New();
	labelStatFilter->SetInput( watersheder->GetOutput() );
	labelStatFilter->Update();

	float meanL = labelStatFilter->GetMean();
	float stdL = labelStatFilter->GetSigma();
	float minL = labelStatFilter->GetMinimum();
	float maxL = labelStatFilter->GetMaximum();

	cout << "Mean : " << meanL << " Std : " << stdL << " Min : " << minL << " Max : " << maxL << endl;

	double* meanArr = new double [ (int)maxL ];
	double* eig0Arr = new double [ (int)maxL ];
	double* eig1Arr = new double [ (int)maxL ];
	double* eig2Arr = new double [ (int)maxL ];

	int* sizeArr = new int[ (int)maxL ];

	for( int i = 0; i < maxL; i++ )
	{
		meanArr[ i ] = 0;
		sizeArr[ i ] = 0;
		eig0Arr[ i ] = 0;
		eig1Arr[ i ] = 0;
		eig2Arr[ i ] = 0;
	}

	for( int x = 0; x < dim3[0]; x++ )
	for( int y = 0; y < dim3[1]; y++ )
	for( int z = 0; z < dim3[2]; z++ )
	{
		ImageType::IndexType idx;
		idx[0] = x;
		idx[1] = y;
		idx[2] = z;

		SegLabelType::IndexType idxSeg;
		idxSeg[0] = x;
		idxSeg[1] = y;
		idxSeg[2] = z;

		unsigned short val = input->GetPixel( idx );
		unsigned long lVal = waterShedLabel->GetPixel( idxSeg );
		double eig0Val = ( eigImg0->GetPixel( idx ) - minEig0 ) * rEig0Inten;
		double eig1Val = ( eigImg1->GetPixel( idx ) - minEig1 ) * rEig1Inten;
		double eig2Val = ( eigImg2->GetPixel( idx ) - minEig2 ) * rEig2Inten;

		int arrIdx = (int)( lVal - 1 );

		double cMean = meanArr[ arrIdx ];
		int cSize = sizeArr[ arrIdx ] ;

		meanArr[ arrIdx ] = cMean + val;
		sizeArr[ arrIdx ] = cSize + 1;
		eig0Arr[ arrIdx ] = eig0Arr[ arrIdx ] + eig0Val;
		eig1Arr[ arrIdx ] = eig1Arr[ arrIdx ] + eig1Val;
		eig2Arr[ arrIdx ] = eig2Arr[ arrIdx ] + eig2Val;
	}

	vector< double > vecMean;
	vecMean.clear();
	vector< double > vecEig0;
	vecEig0.clear();
	vector< double > vecEig1;
	vecEig1.clear();
	vector< double > vecEig2;
	vecEig2.clear();

	vector< int > vecIdx;
	vecIdx.clear();

	int checkSize = 0;
	int check0 = 0;
	int checkTooBig = 0;
	int checkGood = 0;

	for( int i = 0; i < maxL; i++ )
	{
		double meanVal = meanArr[ i ];
		double sizeVal = (double) sizeArr[ i ] ;

		if( sizeVal == 0 )
		{
			check0++;
			continue;
		}
		else if( sizeVal > dim3[0] * dim3[1] * dim3[2] / 3 )
		{
			checkTooBig++;
			continue;
		}
		else
		{
			meanArr[ i ] = meanVal / sizeVal;
			eig0Arr[ i ] = eig0Arr[ i ] / sizeVal;
			eig1Arr[ i ] = eig1Arr[ i ] / sizeVal;
			eig2Arr[ i ] = eig2Arr[ i ] / sizeVal;

			vecMean.push_back( meanArr[ i ] );
			vecIdx.push_back( i );
			vecEig0.push_back( eig0Arr[ i ] );
			vecEig1.push_back( eig1Arr[ i ] );
			vecEig2.push_back( eig2Arr[ i ] );
			checkSize += sizeVal;
			checkGood++;
		}
	}

	if( checkSize == dim3[0] * dim3[1] * dim3[2] )
		cout << "Correct Size -" << "Size 0 : " << check0 << " Size Too big : " << checkTooBig << " Good Size : " << checkGood << endl;
	else
		cout << "Incorrect Size : " << checkSize << "Original Size : " << dim3[0] * dim3[1] * dim3[2] << endl;

	VTK_CREATE( vtkDoubleArray, kMeanArr );
	kMeanArr->SetNumberOfComponents( 1 );
	kMeanArr->SetName( "Mean" );
	kMeanArr->SetNumberOfTuples( vecMean.size() );

	VTK_CREATE( vtkDoubleArray, kEig0Arr );
	kEig0Arr->SetNumberOfComponents( 1 );
	kEig0Arr->SetName( "Eig0" );
	kEig0Arr->SetNumberOfTuples( vecMean.size() );

	VTK_CREATE( vtkDoubleArray, kEig1Arr );
	kEig1Arr->SetNumberOfComponents( 1 );
	kEig1Arr->SetName( "Eig1" );
	kEig1Arr->SetNumberOfTuples( vecMean.size() );

	VTK_CREATE( vtkDoubleArray, kEig2Arr );
	kEig2Arr->SetNumberOfComponents( 1 );
	kEig2Arr->SetName( "Eig2" );
	kEig2Arr->SetNumberOfTuples( vecMean.size() );

	ofstream out2(strRes + "MeanEigSpace.csv");

	for( int i = 0; i < vecMean.size(); i++ )
	{
		kMeanArr->SetValue( i, vecMean.at( i ) );
		kEig0Arr->SetValue( i, vecEig0.at( i ) );
		kEig1Arr->SetValue( i, vecEig1.at( i ) );
		kEig2Arr->SetValue( i, vecEig2.at( i ) );

		out2 << vecIdx.at( i ) << "," << vecMean.at( i ) << ","  << vecEig0.at( i ) << "," << vecEig1.at( i ) << "," << vecEig2.at( i ) << endl;
	}
	out2.close();

	kMeanArr->Modified();
	kEig0Arr->Modified();
	kEig1Arr->Modified();
	kEig2Arr->Modified();

	VTK_CREATE( vtkTable, kMeanInput );
	kMeanInput->AddColumn( kMeanArr );
	kMeanInput->AddColumn( kEig0Arr );
	kMeanInput->AddColumn( kEig1Arr );
	kMeanInput->AddColumn( kEig2Arr );
	kMeanInput->Update();

	VTK_CREATE( vtkKMeansStatistics, kMean );
	kMean->SetInput( vtkStatisticsAlgorithm::INPUT_DATA, kMeanInput );
	kMean->SetColumnStatus( kMeanInput->GetColumnName( 0 ), 1 );
	kMean->SetColumnStatus( kMeanInput->GetColumnName( 1 ), 1 );
	kMean->SetColumnStatus( kMeanInput->GetColumnName( 2 ), 1 );
	kMean->SetColumnStatus( kMeanInput->GetColumnName( 3 ), 1 );
	kMean->RequestSelectedColumns();

//	int nCluster1 = 3;
//
//	for( int i = 1; i < 5; i++)
//	{
//		kMean->SetAssessOption( true );
//		kMean->SetDefaultNumberOfClusters( nCluster1 );
//		kMean->Update();
//
//		vector< vector < double > > meanData;
//		vector< vector < double > > eig0Data;
//		vector< vector < double > > eig1Data;
//		vector< vector < double > > eig2Data;
//
//		for( int i = 0; i < nCluster1; i++ )
//		{
//			meanData.push_back( vector<double> () );
//			eig0Data.push_back( vector<double> () );
//			eig1Data.push_back( vector<double> () );
//			eig2Data.push_back( vector<double> () );
//		}
//
//		double meanIntenCluster[ 3 ] = { 0, 0, 0 }; //, 0 };
//		int cntCluster[ 3 ] = { 0, 0, 0 }; //, 0 };
//
//		double* meanIntenArr = new double [ nCluster1 ];
//		double* meanEig0Arr = new double [ nCluster1 ];
//		double* meanEig1Arr = new double [ nCluster1 ];
//		double* meanEig2Arr = new double [ nCluster1 ];
//
//		for(int i = 0; i < vecMean.size(); i++ )
//		{
//			vtkVariant v = kMean->GetOutput()->GetValue( i, kMean->GetOutput()->GetNumberOfColumns() - 1 );
//
//			int iCluster = v.ToInt();
//
//			meanData[ iCluster ].push_back( vecMean.at( i ) );
//			eig0Data[ iCluster ].push_back( vecEig0.at( i ) );
//			eig1Data[ iCluster ].push_back( vecEig1.at( i ) );
//			eig2Data[ iCluster ].push_back( vecEig2.at( i ) );
//
//			meanIntenCluster[ iCluster ] = meanIntenCluster[ iCluster ] + vecMean.at( i );
//			cntCluster[ iCluster ] = cntCluster[ iCluster ] + 1;
//		}
//	}

	int nCluster1 = 3;
	kMean->SetAssessOption( true );
	kMean->SetDefaultNumberOfClusters( nCluster1 );
	kMean->Update();

	vector < int > vecCluster;
	vecCluster.clear();

	cout << "# of Rows : " << kMean->GetOutput()->GetNumberOfRows() << endl;
	cout << "# of Vectors : " << vecMean.size() << endl;

	ofstream outMean( strRes + "_kmean_Mean.csv" );
	ofstream outEig0( strRes + "_kmean_Eig0.csv" );
	ofstream outEig1( strRes + "_kmean_Eig1.csv" );
	ofstream outEig2( strRes + "_kmean_Eig2.csv" );

	vector< vector < double > > meanData;
	vector< vector < double > > eig0Data;
	vector< vector < double > > eig1Data;
	vector< vector < double > > eig2Data;

	for( int i = 0; i < nCluster1; i++ )
	{
		meanData.push_back( vector<double> () );
		eig0Data.push_back( vector<double> () );
		eig1Data.push_back( vector<double> () );
		eig2Data.push_back( vector<double> () );
	}

	double meanIntenCluster[ 3 ] = { 0, 0, 0 }; //, 0 };
	int cntCluster[ 3 ] = { 0, 0, 0 }; //, 0 };

	for(int i = 0; i < vecMean.size(); i++ )
	{
		vtkVariant v = kMean->GetOutput()->GetValue( i, kMean->GetOutput()->GetNumberOfColumns() - 1 );
		int iCluster = v.ToInt();
		vecCluster.push_back( iCluster );

		meanData[ iCluster ].push_back( vecMean.at( i ) );
		eig0Data[ iCluster ].push_back( vecEig0.at( i ) );
		eig1Data[ iCluster ].push_back( vecEig1.at( i ) );
		eig2Data[ iCluster ].push_back( vecEig2.at( i ) );

		meanIntenCluster[ iCluster ] = meanIntenCluster[ iCluster ] + vecMean.at( i );
		cntCluster[ iCluster ] = cntCluster[ iCluster ] + 1;
	}


	for( int i = 0; i < nCluster1; i++ )
	{
		for( int j = 0; j < meanData[ i ].size(); j++ )
		{
			outMean << meanData[ i ].at( j ) << ",";
			outEig0 << eig0Data[ i ].at( j ) << ",";
			outEig1 << eig1Data[ i ].at( j ) << ",";
			outEig2 << eig2Data[ i ].at( j ) << ",";
		}
		outMean << endl;
		outEig0 << endl;
		outEig1 << endl;
		outEig2 << endl;

	}

	outMean.close();
	outEig0.close();
	outEig1.close();
	outEig2.close();


	for( int i = 0; i < nCluster1; i++ )
	{
		meanIntenCluster[ i ] = meanIntenCluster[ i ] / (double) cntCluster[ i ];
	}

	// Find Cluster with minimum intensity
	int idxMinCl = 0;
	double minMeanCl = 10000000;

	for( int i = 0; i < nCluster1; i++ )
	{
		if( minMeanCl > meanIntenCluster[ i ] )
		{
			minMeanCl = meanIntenCluster[ i ];
			idxMinCl = i;
		}
	}

	LabelType::Pointer kMeanLabel = LabelType::New();
	kMeanLabel->SetRegions( input->GetLargestPossibleRegion() );
	kMeanLabel->SetSpacing( input->GetSpacing() );
	kMeanLabel->SetOrigin( input->GetOrigin() );
	kMeanLabel->Allocate();

	for( int x = 0; x < dim3[0]; x++ )
	for( int y = 0; y < dim3[1]; y++ )
	for( int z = 0; z < dim3[2]; z++ )
	{
		LabelType::IndexType idx;
		idx[0] = x;
		idx[1] = y;
		idx[2] = z;

		unsigned long lVal = waterShedLabel->GetPixel( idx );

		int arrIdx = (int)( lVal - 1 );

		kMeanLabel->SetPixel( idx, 0 );

		for( int j = 0; j < vecIdx.size(); j++ )
		{
			if( arrIdx == vecIdx.at( j ) )
			{
//				kMeanLabel->SetPixel( idx, vecCluster.at( j ) );
				if( vecCluster.at( j ) == idxMinCl )
				{
					kMeanLabel->SetPixel( idx, 2 );
					break;
				}
			}
		}
	}

	kMeanLabel->Update();

	labelWriter->SetInput( kMeanLabel );
	labelWriter->SetFileName(strRes + "kMeanLabel.nrrd" );
	labelWriter->Update();

	CCAFilterType::Pointer ccaFilter = CCAFilterType::New();
	ccaFilter->SetInput( kMeanLabel );
	ccaFilter->Update();

	RelabelFilterType::Pointer sizeFilter = RelabelFilterType::New();
	RelabelFilterType::ObjectSizeType  minSize = 0.4;
	sizeFilter->SetMinimumObjectSize( minSize );
	sizeFilter->SetInput( ccaFilter->GetOutput() );
	sizeFilter->Update();

	labelWriter->SetInput( sizeFilter->GetOutput() );
	labelWriter->SetFileName( strRes + "kMeanCCA.nrrd" );
	labelWriter->Update();




	// 2nd Clustering Begin
	vector< unsigned short > vecBlobInten;
	vecBlobInten.clear();

	vector< ImageType::IndexType > vecBlobIdx;
	vecBlobIdx.clear();

	vector< double > vecBlobEig0;
	vecBlobEig0.clear();

	vector< double > vecBlobEig1;
	vecBlobEig1.clear();

	vector< double > vecBlobEig2;
	vecBlobEig2.clear();

	for( int x = 0; x < dim3[0]; x++ )
	for( int y = 0; y < dim3[1]; y++ )
	for( int z = 0; z < dim3[2]; z++ )
	{
		ImageType::IndexType idx;
		idx[0] = x;
		idx[1] = y;
		idx[2] = z;

		unsigned char kL = kMeanLabel->GetPixel( idx );

		if( kL == 2 )
		{
			vecBlobInten.push_back( input->GetPixel( idx ) );
			vecBlobEig0.push_back( abs( eigImg0->GetPixel( idx ) ) * rEig0Inten );
			vecBlobEig1.push_back( abs( eigImg1->GetPixel( idx ) ) * rEig1Inten );
			vecBlobEig2.push_back( abs( eigImg2->GetPixel( idx ) ) * rEig2Inten );
			vecBlobIdx.push_back( idx );
		}
	}

	VTK_CREATE( vtkDoubleArray, blobIntenArray );
	blobIntenArray->SetNumberOfComponents( 1 );
	blobIntenArray->SetName( "Intensity" );
	blobIntenArray->SetNumberOfTuples( vecBlobInten.size() );

	VTK_CREATE( vtkDoubleArray, blobBlobnessArray );
	blobBlobnessArray->SetNumberOfComponents( 1 );
	blobBlobnessArray->SetName( "Blobness" );
	blobBlobnessArray->SetNumberOfTuples( vecBlobInten.size() );

	VTK_CREATE( vtkDoubleArray, blobVesselnessArray );
	blobVesselnessArray->SetNumberOfComponents( 1 );
	blobVesselnessArray->SetName( "Vesselness" );
	blobVesselnessArray->SetNumberOfTuples( vecBlobInten.size() );

	for( int i = 0; i < vecBlobInten.size(); i++ )
	{
		double eigVal0 = vecBlobEig0.at( i );
		double eigVal1 = vecBlobEig1.at( i );
		double eigVal2 = vecBlobEig2.at( i );
		RealImageType::IndexType idx = vecBlobIdx.at( i );


		double blobness = blobImg->GetPixel( idx ) * rBlobInten;
		double vesselness = vesselImg->GetPixel( idx ) * rVesInten;

		blobIntenArray->SetValue( i, vecBlobInten.at( i ) );
		blobVesselnessArray->SetValue( i, vesselness );
		blobBlobnessArray->SetValue( i, blobness );
	}

	blobIntenArray->Modified();
	blobBlobnessArray->Modified();
	blobVesselnessArray->Modified();

	VTK_CREATE( vtkTable, tableBlobInten );
	tableBlobInten->AddColumn( blobIntenArray );
	tableBlobInten->AddColumn( blobBlobnessArray );
	tableBlobInten->AddColumn( blobVesselnessArray );
	tableBlobInten->Update();

	// 2nd KMeans
	VTK_CREATE( vtkKMeansStatistics, kMean2 );
	kMean2->SetInput( vtkStatisticsAlgorithm::INPUT_DATA, tableBlobInten );
	kMean2->SetColumnStatus( tableBlobInten->GetColumnName( 0 ), 1 );
	kMean2->SetColumnStatus( tableBlobInten->GetColumnName( 1 ), 1 );
	kMean2->SetColumnStatus( tableBlobInten->GetColumnName( 2 ), 1 );
//	kMean2->SetColumnStatus( tableBlobInten->GetColumnName( 3 ), 1 );

	kMean2->RequestSelectedColumns();
	kMean2->SetAssessOption( true );
	kMean2->SetDefaultNumberOfClusters( 2 );
	kMean2->Update();
//	kMean->GetOutput()->Dump();

	vector < int > vecCluster2;
	vecCluster2.clear();

	LabelType::Pointer kMeanLabel2 = LabelType::New();
	kMeanLabel2->SetRegions( input->GetLargestPossibleRegion() );
	kMeanLabel2->SetSpacing( input->GetSpacing() );
	kMeanLabel2->SetOrigin( input->GetOrigin() );
	kMeanLabel2->Allocate();

	double meanI2ndCl0;
	double meanI2ndCl1;
	double cnt2ndCl0;
	double cnt2ndCl1;

	for( int x = 0; x < dim3[0]; x++ )
	for( int y = 0; y < dim3[1]; y++ )
	for( int z = 0; z < dim3[2]; z++ )
	{
		ImageType::IndexType idx;
		idx[0] = x;
		idx[1] = y;
		idx[2] = z;

		kMeanLabel2->SetPixel( idx, 0 );
	}

	for( int i = 0; i < vecBlobInten.size(); i++ )
	{
		ImageType::IndexType idx = vecBlobIdx.at( i );
		vtkVariant v = kMean2->GetOutput()->GetValue( i, kMean2->GetOutput()->GetNumberOfColumns() - 1 );
		int iCluster = v.ToInt() + 1;

		kMeanLabel2->SetPixel( idx, iCluster );

		if( iCluster == 1 )
		{
			meanI2ndCl0 += vecBlobInten.at( i );
			cnt2ndCl0++;
		}
		else
		{
			meanI2ndCl1 += vecBlobInten.at( i );
			cnt2ndCl1++;
		}
	}

	kMeanLabel2->Update();

	labelWriter->SetInput( kMeanLabel2 );
	labelWriter->SetFileName( strRes + "kMeanLabel2nd.nrrd" );
	labelWriter->Update();

	// Write Only Blob Candidates
	meanI2ndCl0 /= cnt2ndCl0;
	meanI2ndCl1 /= cnt2ndCl1;

	int min2ndIdx = 0;

	if( meanI2ndCl0 > meanI2ndCl1 )
		min2ndIdx = 1;
	else
		min2ndIdx = 2;

	LabelType::Pointer kMeanLabel3 = LabelType::New();
	kMeanLabel3->SetRegions( input->GetLargestPossibleRegion() );
	kMeanLabel3->SetSpacing( input->GetSpacing() );
	kMeanLabel3->SetOrigin( input->GetOrigin() );
	kMeanLabel3->Allocate();

	for( int x = 0; x < dim3[0]; x++ )
	for( int y = 0; y < dim3[1]; y++ )
	for( int z = 0; z < dim3[2]; z++ )
	{
		ImageType::IndexType idx;
		idx[0] = x;
		idx[1] = y;
		idx[2] = z;

		kMeanLabel3->SetPixel( idx, 0 );
	}

	for( int i = 0; i < vecBlobInten.size(); i++ )
	{
		ImageType::IndexType idx = vecBlobIdx.at( i );
		vtkVariant v = kMean2->GetOutput()->GetValue( i, kMean2->GetOutput()->GetNumberOfColumns() - 1 );
		int iCluster = v.ToInt() ;

		if( iCluster == min2ndIdx )
			kMeanLabel3->SetPixel( idx, iCluster );
	}

	kMeanLabel3->Update();

	labelWriter->SetInput( kMeanLabel3 );
	labelWriter->SetFileName( strRes + "kMeanLabel2ndBlobCandidates.nrrd" );
	labelWriter->Update();

	CCAFilterType::Pointer ccaFilter2nd = CCAFilterType::New();
	ccaFilter2nd->SetInput( kMeanLabel3 );
	ccaFilter2nd->Update();

	SegRelabelFilterType::Pointer sizeFilter2nd = SegRelabelFilterType::New();
	sizeFilter2nd->SetMinimumObjectSize( 0.4 );
	sizeFilter2nd->SetInput( ccaFilter2nd->GetOutput() );
	sizeFilter2nd->Update();

	SegLabelType::Pointer blobCCALabel = sizeFilter2nd->GetOutput();

	segLabelWriter->SetInput( blobCCALabel );
	segLabelWriter->SetFileName( strRes + "kMeanLabel2ndBlobCandidatesCCA.nrrd" );
	segLabelWriter->Update();

//	// 2nd Cluster Result Shape Refinement
//	// Get Cubic OBB
//	int nObj = sizeFilter2nd->GetNumberOfObjects();
//
//	cout << "# of Objects : " << nObj << endl;
//
//	int* lx = new int[ nObj ];
//	int* hx = new int[ nObj ];
//
//	int* ly = new int[ nObj ];
//	int* hy = new int[ nObj ];
//
//	int* lz = new int[ nObj ];
//	int* hz = new int[ nObj ];
//
//	for( int i = 0; i < nObj; i++ )
//	{
//		lx[ i ] = 100000;
//		ly[ i ] = 100000;
//		lz[ i ] = 100000;
//		hx[ i ] = 0;
//		hy[ i ] = 0;
//		hz[ i ] = 0;
//	}
//
//	for( int x = 0; x < dim3[0]; x++ )
//	for( int y = 0; y < dim3[1]; y++ )
//	for( int z = 0; z < dim3[2]; z++ )
//	{
//		SegLabelType::IndexType idx;
//		idx[0] = x;
//		idx[1] = y;
//		idx[2] = z;
//
//		unsigned long segLabelVal = blobCCALabel->GetPixel( idx );
//
//		if( segLabelVal == 0 )
//			continue;
//
//		int labIdx = segLabelVal - 1;
//
//		if( labIdx < 0 )
//			continue;
//
//		if( lx[ labIdx ] > x ) lx[ labIdx ] = x;
//		if( ly[ labIdx ] > y ) ly[ labIdx ] = y;
//		if( lz[ labIdx ] > z ) lz[ labIdx ] = z;
//
//		if( hx[ labIdx ] < x ) hx[ labIdx ] = x;
//		if( hy[ labIdx ] < y ) hy[ labIdx ] = y;
//		if( hz[ labIdx ] < z ) hz[ labIdx ] = z;
//	}
//
//	// Calculate the ratio between width and height and the ratio of object volume in OBB
//	vector< int > vecFlaws;
//	vecFlaws.clear();
//
//	for( int i = 0; i < nObj; i++ )
//	{
//		double ratioMax = 0;
//		double ratioXY = (double)( hx[i] - lx[i] + 1) / (double)( hy[i] - ly[i] + 1);
//		double ratioXZ = (double)( hx[i] - lx[i] + 1) / (double)( hz[i] - lz[i] + 1);
//		double ratioYZ = (double)( hy[i] - ly[i] + 1) / (double)( hz[i] - lz[i] + 1);
//		double ratioYX = 1.0 / ratioXY;
//		double ratioZX = 1.0 / ratioXZ;
//		double ratioZY = 1.0 / ratioYZ;
//
//		double whR[6] = { ratioXY, ratioXZ, ratioYZ, ratioYX, ratioZX, ratioZY };
//
//		for( int j = 0; j < 6; j++ )
//		{
//			if( ratioMax < whR[ j ] ) ratioMax = whR[ j ];
//		}
//
//
//		double obbVol = ( hx[ i ] - lx[ i ] + 1 ) * ( hy[ i ] - ly[ i ] + 1 ) * ( hz[ i ] - lz[ i ] + 1 );
//		if( obbVol <= 0 ) obbVol = 1;
//
//		double ratioVol = (double) sizeFilter2nd->GetSizeOfObjectsInPixels()[ i ] / obbVol;
//
//		cout << "Index : " << i << " Volume Ratio : " << ratioVol << " Line Ratio : " << ratioMax << endl;
//
//		if( ratioMax > 4.0 || ratioVol < 0.1 )
//		{
//			vecFlaws.push_back( i );
//		}
//	}
//
//	SegLabelType::Pointer kMeanLabelRefined = SegLabelType::New();
//	kMeanLabelRefined->SetRegions( input->GetLargestPossibleRegion() );
//	kMeanLabelRefined->SetSpacing( input->GetSpacing() );
//	kMeanLabelRefined->SetOrigin( input->GetOrigin() );
//	kMeanLabelRefined->Allocate();
//
//	for( int x = 0; x < dim3[0]; x++ )
//	for( int y = 0; y < dim3[1]; y++ )
//	for( int z = 0; z < dim3[2]; z++ )
//	{
//		ImageType::IndexType idx;
//		idx[0] = x;
//		idx[1] = y;
//		idx[2] = z;
//
//		unsigned long segLabelVal = blobCCALabel->GetPixel( idx );
//		kMeanLabelRefined->SetPixel( idx, segLabelVal );
//
//		for( int i = 0; i < vecFlaws.size(); i++ )
//		{
//			if( segLabelVal == vecFlaws.at( i ) )
//			{
//				kMeanLabelRefined->SetPixel( idx, 0 );
//				break;
//			}
//		}
//	}
//
//	kMeanLabelRefined->Update();
//
//	segLabelWriter->SetInput( kMeanLabelRefined );
//	segLabelWriter->SetFileName( strRes + "kMeanLabel2ndBlobRefined.nrrd" );
//	segLabelWriter->Update();
//
//	// Termination
//	delete [] lx;
//	delete [] ly;
//	delete [] lz;
//	delete [] hx;
//	delete [] hy;
//	delete [] hz;


//	/***********************************************
//	 * Blobs are a problem in the image. They split
//	 * fibers into pieces. Most of the blobs have
//	 * the same label as fibers, but are isolated.
//	 * If we can segment such blobs, we can dilate
//	 * them and only them so they touch the rest
//	 * of the skeleton they might have split.
//	 * If I is the binary image and B the binary
//	 * image with only the blobs, D the
//	 * dilatation operator, and K the skeletonization
//	 * operator, K( I or D(B) ) should get rid of
//	 * the split artifact.
//	 ***********************************
//	***********/
//
//
//	std::cout << " ***** Skeletonization " << std::endl;
//	// We give the binary output of the last filter
//
//	BinaryThinningFilterType::Pointer thinningFilter = BinaryThinningFilterType::New();
//	thinningFilter->SetInput( lineShapeFilter->GetOutput() );
//	thinningFilter->Update();
//
//	writer->SetInput( thinningFilter->GetOutput() );
//	writer->SetFileName( str + "_skeleton.nrrd" );
//	writer->Update();
//
//	std::cout << " ***** Vesselness Skeleton " << std::endl;
//	// We get an image of the skeleton masked with the vesselness image
//
//	MultiplyFilterType::Pointer multiplyFilter = MultiplyFilterType::New();
//	multiplyFilter->SetInput1( thinningFilter->GetOutput() );
//	multiplyFilter->SetInput2( lineShapeFilter->GetVesselnessOutput() );
//	multiplyFilter->Update();
//
//	realWriter->SetInput( multiplyFilter->GetOutput() );
//	realWriter->SetFileName( str + "_skeleton_vesselness.nrrd" );
//	realWriter->Update();
//
//	std::cout << " ***** Skeleton To Graph Filter " << std::endl;
//	// From the skeleton, we build a graph
//
//	SkeletonToGraphFilterType::Pointer skeletonToGraphFilter = SkeletonToGraphFilterType::New();
//	skeletonToGraphFilter->SetInput( thinningFilter->GetOutput() );
//	skeletonToGraphFilter->SetVesselnessInput( lineShapeFilter->GetVesselnessOutput() );
//	skeletonToGraphFilter->Update();
//
//	// Identifier output of the skeleton to graph
//	// This is a skeleton on which values are the number of the nodes
//	itk::ImageFileWriter< SkeletonToGraphFilterType::IdentifierImageType >::Pointer iwriter = itk::ImageFileWriter< SkeletonToGraphFilterType::IdentifierImageType >::New();
//	iwriter->SetInput( skeletonToGraphFilter->GetIdentifierOutput() );
//	iwriter->SetFileName( str + "_identifier.nrrd" );
//	iwriter->Update();
//
//	std::cout << " ***** Fiber extraction Graph Filter  " << std::endl;
//	// It takes the graph as input and processes it
//	// This version works with straight fibers
//
//	FiberExtractionFilterType::Pointer fiberExtractionFilter = FiberExtractionFilterType::New();
//	fiberExtractionFilter->SetInput( skeletonToGraphFilter->GetOutput() );
//	fiberExtractionFilter->Update();
//
//	std::cout << " ***** Converting graph to image " << std::endl;
//	// To see what we get
//	// If you want to see any values, go into the GraphToSkeletonFilter
//	// file, and change what you write (very small class).
//
//	GraphToSkeletonFilterType::Pointer graphToSkeletonFilter = GraphToSkeletonFilterType::New();
//	graphToSkeletonFilter->SetInput( fiberExtractionFilter->GetOutput() );
//	graphToSkeletonFilter->SetRegion( reader->GetOutput()->GetLargestPossibleRegion() );
//	graphToSkeletonFilter->Update();
//
//	// Writing what we want to see
//	realWriter->SetInput( graphToSkeletonFilter->GetOutput() );
//	realWriter->SetFileName( str + "_fiber_extraction.nrrd" );
//	realWriter->Update();

	cout << "Process Done" << endl;
	return EXIT_SUCCESS;
}
