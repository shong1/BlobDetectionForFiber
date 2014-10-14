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

typedef itk::ImageFileWriter< RealImageType >							RealWriterType;

typedef itk::SkeletonToGraphFilter< ImageType, GraphType, RealImageType, ImageType >	SkeletonToGraphFilterType;
typedef itk::GraphToSkeletonImageFilter< GraphType, RealImageType >			GraphToSkeletonFilterType;
typedef itk::FiberExtractionGraphFilter< GraphType >					FiberExtractionFilterType;

typedef itk::MultiplyImageFilter< ImageType, RealImageType, RealImageType > MultiplyFilterType;

typedef itk::ImageDuplicator< RealImageType >	DuplicatorType;

typedef itk::MinimumMaximumImageCalculator< ImageType >					MinMaxFilterType;
typedef itk::StatisticsImageFilter< ImageType > 						StaticFilterType;
typedef itk::ScalarImageKmeansImageFilter< ImageType, LabelType > 		KMeansType;
typedef itk::StatisticsImageFilter< SegLabelType > 						LabelStaticFilterType;

typedef itk::ConnectedComponentImageFilter< LabelType, ImageType >		CCAFilterType;
typedef itk::RelabelComponentImageFilter< ImageType, LabelType > 		RelabelFilterType;
typedef itk::BinaryThresholdImageFilter< ImageType, LabelType > 		BinaryFilterType;
typedef itk::BinaryBallStructuringElement< LabelType::PixelType, 3 >	StructureElementType;
typedef itk::BinaryErodeImageFilter< LabelType, LabelType, StructureElementType > 			BinErodeFilterType;
typedef itk::GrayscaleDilateImageFilter< ImageType, ImageType, StructureElementType >		GrayDilateFilterType;
typedef itk::SubtractImageFilter< ImageType, ImageType > 				SubtractorType;

typedef itk::GradientMagnitudeImageFilter< ImageType, RealImageType > 	GradientMagnitudeFilterType;
typedef itk::WatershedImageFilter< RealImageType > 						WatershedFilterType;

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
	std::string str = "/home/shong/Research/FibersInConcrete/Data/StiffFibers_cropped";

	WriterType::Pointer writer = WriterType::New();
	RealWriterType::Pointer realWriter = RealWriterType::New();
	LabelWriterType::Pointer labelWriter = LabelWriterType::New();
	SegLabelWriterType::Pointer segLabelWriter = SegLabelWriterType::New();

	/*
	 *	Reading input image
	 */

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

	float statMin = meanI - 3 * stdI;
	float statMax = meanI + 3 * stdI;
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


		unsigned short val = input->GetPixel( idx );

		if( val < statMin )
			input->SetPixel( idx, statMin );
		else if( val > statMax )
			input->SetPixel( idx, statMax );
	}

	input->Update();

	writer->SetInput( input );
	writer->SetFileName( str + "_IntensityCut.nrrd" );
	writer->Update();


	std::cout << " ***** LineShapeFilter " << std::endl;

	// This filter does a hessian + gaussian filter
	// We have then access to eigen values
	LineShapeFilterType::Pointer lineShapeFilter = LineShapeFilterType::New();

	lineShapeFilter->SetInput( reader->GetOutput() );
	// this parameter has to be about the width of a fiber
	lineShapeFilter->SetSigma( 2.0f );
	lineShapeFilter->SetExtractBrightLine( false );
	lineShapeFilter->EigenValuesExtractionOn();
	lineShapeFilter->LabelImageOn();
	lineShapeFilter->Update();
	// Binary image of the fiber-like shapes
	writer->SetInput( lineShapeFilter->GetBinaryOutput() );
	writer->SetFileName( str + "_binary.nrrd" );
	writer->Update();

	// Label image depending of the sign of the 3 sorted eigenvalues
	writer->SetInput( lineShapeFilter->GetLabelOutput() );
	writer->SetFileName( str + "_label.nrrd" );
	writer->Update();

	RealImageType::Pointer eigImg0 = lineShapeFilter->GetEigenValuesOutput( 0 );
	RealImageType::Pointer eigImg1 = lineShapeFilter->GetEigenValuesOutput( 1 );
	RealImageType::Pointer eigImg2 = lineShapeFilter->GetEigenValuesOutput( 2 );

	// Watershed Segmentation for Blob Detection
	GradientMagnitudeFilterType::Pointer gradientFilter = GradientMagnitudeFilterType::New();
	gradientFilter->SetInput( input );
	gradientFilter->Update();

//	realWriter->SetInput( gradientFilter->GetOutput() );
//	realWriter->SetFileName( str + "_gradMag.nrrd" );
//	realWriter->Update();

	WatershedFilterType::Pointer watersheder = WatershedFilterType::New();
	watersheder->SetInput( gradientFilter->GetOutput() );
	watersheder->SetThreshold( 0.001 );
	watersheder->SetLevel( 0.2 );
	watersheder->Update();

	segLabelWriter->SetInput( watersheder->GetOutput() );
	segLabelWriter->SetFileName( str + "_watershed_001_2.nrrd" );
	segLabelWriter->Update();

	LabelStaticFilterType::Pointer labelStatFilter = LabelStaticFilterType::New();
	labelStatFilter->SetInput( watersheder->GetOutput() );
	labelStatFilter->Update();

	SegLabelType::Pointer waterShedLabel = watersheder->GetOutput();

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

	for( int i = 1; i < maxL; i++ )
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
		double eig0Val = eigImg0->GetPixel( idx );
		double eig1Val = eigImg1->GetPixel( idx );
		double eig2Val = eigImg2->GetPixel( idx );


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

	for( int i = 0; i < maxL; i++ )
	{
		double meanVal = meanArr[ i ];
		double sizeVal = (double) sizeArr[ i ] ;

		if( sizeVal == 0 )
			continue;
		else if( sizeVal > dim3[0] * dim3[1] * dim3[2] / 3 )
			continue;
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
		}
	}

	cout << "Vector Size : " << vecMean.size() << endl;

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

	ofstream out2(str + "Mean.csv");

	for( int i = 0; i < vecMean.size(); i++ )
	{
		kMeanArr->SetValue( i, vecMean.at( i ) );
		kEig0Arr->SetValue( i, vecEig0.at( i ) );
		kEig1Arr->SetValue( i, vecEig1.at( i ) );
		kEig2Arr->SetValue( i, vecEig2.at( i ) );

		out2 << vecMean.at( i ) << "," << vecIdx.at( i ) << endl;
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
	kMean->SetAssessOption( true );
	kMean->SetDefaultNumberOfClusters( 4 );
	kMean->Update();
//	kMean->GetOutput()->Dump();

	vector < int > vecCluster;
	vecCluster.clear();

	cout << "# of Rows : " << kMean->GetOutput()->GetNumberOfRows() << endl;
	cout << "# of Vectors : " << vecMean.size() << endl;

	ofstream out( str + "kmean.csv" );

	double meanIntenCluster[ 4 ] = { 0, 0, 0, 0 };
	int cntCluster[ 4 ] = { 0, 0, 0, 0 };

	for(int i = 0; i < vecMean.size(); i++ )
	{
		vtkVariant v = kMean->GetOutput()->GetValue( i, kMean->GetOutput()->GetNumberOfColumns() - 1 );
		int iCluster = v.ToInt();
		vecCluster.push_back( iCluster );
		out << iCluster;
		out << ",";
		out << vecMean.at(i);
		out << endl;

		meanIntenCluster[ iCluster ] = meanIntenCluster[ iCluster ] + vecMean.at( i );
		cntCluster[ iCluster ] = cntCluster[ iCluster ] + 1;
	}

	for( int i = 0; i < 4; i++ )
	{
		meanIntenCluster[ i ] = meanIntenCluster[ i ] / (double) cntCluster[ i ];
	}

	// Find Cluster with minimum intensity
	int idxMinCl = 0;
	double minMeanCl = 10000000;

	for( int i = 0; i < 4; i++ )
	{
		if( minMeanCl > meanIntenCluster[ i ] )
		{
			minMeanCl = meanIntenCluster[ i ];
			idxMinCl = i;
		}
	}
	out.close();

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
	labelWriter->SetFileName(str + "kMeanLabel.nrrd" );
	labelWriter->Update();

	CCAFilterType::Pointer ccaFilter = CCAFilterType::New();
	ccaFilter->SetInput( kMeanLabel );
	ccaFilter->Update();

	RelabelFilterType::Pointer sizeFilter = RelabelFilterType::New();
	RelabelFilterType::ObjectSizeType  minSize = 5;
	sizeFilter->SetMinimumObjectSize( minSize );
	sizeFilter->SetInput( ccaFilter->GetOutput() );
	sizeFilter->Update();

	labelWriter->SetInput( sizeFilter->GetOutput() );
	labelWriter->SetFileName( str + "kMeanCCA.nrrd" );
	labelWriter->Update();

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
			vecBlobEig0.push_back( eigImg0->GetPixel( idx ) );
			vecBlobEig1.push_back( eigImg1->GetPixel( idx ) );
			vecBlobEig2.push_back( eigImg2->GetPixel( idx ) );
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

	VTK_CREATE( vtkDoubleArray, blobStructurenessArray );
	blobStructurenessArray->SetNumberOfComponents( 1 );
	blobStructurenessArray->SetName( "Structureness" );
	blobStructurenessArray->SetNumberOfTuples( vecBlobInten.size() );

	VTK_CREATE( vtkDoubleArray, blobVesselnessArray );
	blobVesselnessArray->SetNumberOfComponents( 1 );
	blobVesselnessArray->SetName( "Vesselness" );
	blobVesselnessArray->SetNumberOfTuples( vecBlobInten.size() );


	for( int i = 0; i < vecBlobInten.size(); i++ )
	{
		double eigVal0 = vecBlobEig0.at( i );
		double eigVal1 = vecBlobEig1.at( i );
		double eigVal2 = vecBlobEig2.at( i );

		double blobness = abs( eigVal0 ) / sqrt( abs( eigVal1 * eigVal2 ) ) * 20000;
		double structureness = eigVal0 * eigVal0 + eigVal1 * eigVal1 + eigVal2 * eigVal2;
		double vesselness = abs( eigVal1 / eigVal2 ) * 20000;

		blobIntenArray->SetValue( i, vecBlobInten.at( i ) );
		blobBlobnessArray->SetValue( i, blobness );
		blobStructurenessArray->SetValue( i, structureness );
		blobVesselnessArray->SetValue( i, vesselness );
	}

	blobIntenArray->Modified();
	blobBlobnessArray->Modified();
	blobStructurenessArray->Modified();
	blobVesselnessArray->Modified();

	VTK_CREATE( vtkTable, tableBlobInten );
	tableBlobInten->AddColumn( blobIntenArray );
	tableBlobInten->AddColumn( blobBlobnessArray );
	tableBlobInten->AddColumn( blobStructurenessArray );
	tableBlobInten->AddColumn( blobVesselnessArray );
	tableBlobInten->Update();

	// 2nd KMeans
	VTK_CREATE( vtkKMeansStatistics, kMean2 );
	kMean2->SetInput( vtkStatisticsAlgorithm::INPUT_DATA, tableBlobInten );
	kMean2->SetColumnStatus( tableBlobInten->GetColumnName( 0 ), 1 );

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
	labelWriter->SetFileName(str + "kMeanLabel2nd.nrrd" );
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
	labelWriter->SetFileName(str + "kMeanLabel2ndBlobCandidates.nrrd" );
	labelWriter->Update();

	CCAFilterType::Pointer ccaFilter2nd = CCAFilterType::New();
	ccaFilter2nd->SetInput( kMeanLabel3 );
	ccaFilter2nd->Update();

	RelabelFilterType::Pointer sizeFilter2nd = RelabelFilterType::New();
	sizeFilter2nd->SetMinimumObjectSize( minSize );
	sizeFilter2nd->SetInput( ccaFilter2nd->GetOutput() );
	sizeFilter2nd->Update();

	labelWriter->SetInput( sizeFilter2nd->GetOutput() );
	labelWriter->SetFileName( str + "kMeanLabel2ndBlobCandidatesCCA.nrrd" );
	labelWriter->Update();

	int nObj = sizeFilter2nd->GetNumberOfObjects();

	// Find Cluster with minimum intensity
	cout << "Process Done" << endl;


/*
	RealImageType::Pointer eigImg0 = lineShapeFilter->GetEigenValuesOutput( 0 );
	RealImageType::Pointer eigImg1 = lineShapeFilter->GetEigenValuesOutput( 1 );
	RealImageType::Pointer eigImg2 = lineShapeFilter->GetEigenValuesOutput( 2 );

	RealImageType::RegionType region = eigImg0->GetLargestPossibleRegion();
	RealImageType::SizeType dim = region.GetSize();

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

	RealImageType::Pointer vesselImg = RealImageType::New();
	vesselImg->SetRegions( eigImg0->GetLargestPossibleRegion() );
	vesselImg->SetSpacing( eigImg0->GetSpacing() );
	vesselImg->SetOrigin( eigImg0->GetOrigin() );
	vesselImg->Allocate();


	cout << dim[0] << ", " << dim[1] << ", " << dim[2] << endl;
	float alpha = 0.5;
	float beta = 0.5;
	float c = 0.5;

	for( int x = 0; x < dim[0]; x++ )
	for( int y = 0; y < dim[1]; y++ )
	for( int z = 0; z < dim[2]; z++ )
	{
		RealImageType::IndexType idx;
		idx[0] = x;
		idx[1] = y;
		idx[2] = z;
		float lambda1 = eigImg0->GetPixel( idx );
		float lambda2 = eigImg1->GetPixel( idx );
		float lambda3 = eigImg2->GetPixel( idx );

		// eigImg2 > eigImg1 > eigImg0
		// lambda3 > lambda2 > lambda1
		float blVal = abs( lambda1 ) / sqrt( abs( lambda2 * lambda3 ) );
		float sVal = sqrt( lambda1 * lambda1 + lambda2 * lambda2 + lambda3 * lambda3 );
		float aVal = abs( lambda2 / lambda3 );

		float V0 = exp( -( aVal*aVal ) / ( 2 * alpha * alpha ) ) * ( 1 - exp( -(blVal * blVal) / ( 2 * beta * beta ) ) ) * ( 1 - exp( -( sVal * sVal ) / ( 2 * c * c ) ) );

		float Ves = ( 1 - exp( -( aVal*aVal ) / ( 2 * alpha * alpha ) ) ) * ( exp( -(blVal * blVal) / ( 2 * beta * beta ) ) ) * ( 1 - exp( -( sVal * sVal ) / ( 2 * c * c ) ) );
//		float V0 = 1 - exp( ( -blVal * blVal ) / ( 2 * beta * beta ) );
		blobImg->SetPixel( idx, V0 );
		vesselImg->SetPixel( idx, blVal );

		if( V0 > 0.2 && lambda1 < 0 && lambda2 < 0 && lambda3 < 0 )
			blobBin->SetPixel( idx, 255 );
		else
			blobBin->SetPixel( idx, 0 );
	}

	blobImg->Update();
	blobBin->Update();
	vesselImg->Update();

	cout << "Blob Detected?" << endl;

	realWriter->SetInput( blobImg );
	realWriter->SetFileName( str + "blobness.nrrd" );
	realWriter->Update();

	realWriter->SetInput( vesselImg );
	realWriter->SetFileName( str + "vesselnessF.nrrd" );
	realWriter->Update();


	// Label image depending of the sign of the 3 sorted eigenvalues
//	writer->SetInput( blobBin );
//	writer->SetFileName( str + "_label.nrrd" );
//	writer->Update();

	cout << "Blob Written" << endl;

	CCAFilterType::Pointer ccaFilter = CCAFilterType::New();
	ccaFilter->SetInput( blobBin );
	ccaFilter->Update();

	cout << "# Objects : " << ccaFilter->GetObjectCount() << endl;

	RelabelFilterType::Pointer sizeFilter = RelabelFilterType::New();
	RelabelFilterType::ObjectSizeType  minSize = 10;
	sizeFilter->SetInput( ccaFilter->GetOutput () );
	sizeFilter->SetMinimumObjectSize( minSize );
	sizeFilter->Update();

	labelWriter->SetInput( sizeFilter->GetOutput() );
	labelWriter->SetFileName( str + "labelS.nrrd" );
	labelWriter->Update();

	// Vesselness output (fiber-shape likelihood)
	realWriter->SetInput( lineShapeFilter->GetVesselnessOutput() );
	realWriter->SetFileName( str + "_vesselness.nrrd" );
	realWriter->Update();

	realWriter->SetInput( blobImg );
	realWriter->SetFileName( str + "_blobness.nrrd" );
	realWriter->Update();


	// Label image depending of the sign of the 3 sorted eigenvalues
	writer->SetInput( lineShapeFilter->GetLabelOutput() );
	writer->SetFileName( str + "_label.nrrd" );
	writer->Update();

	cout << "Process Done"<< endl;
*/


//
//	KMeansType::Pointer KmeanFilter = KMeansType::New();
//	KmeanFilter->SetInput( input );
//	KmeanFilter->SetUseNonContiguousLabels( true );
//	KmeanFilter->AddClassWithInitialMean( statMin );
//	KmeanFilter->AddClassWithInitialMean( meanI );
//	KmeanFilter->AddClassWithInitialMean( statMax );
//	KmeanFilter->Update();
//
//	KMeansType::ParametersType estimatedMean = KmeanFilter->GetFinalMeans();
//
//	int nClass = estimatedMean.Size();
//
//	cout << "# of clusters : " << nClass << endl;
//
//	labelWriter->SetInput( KmeanFilter->GetOutput() );
//	labelWriter->SetFileName( str + "kmeansLabel.nrrd" );
//	labelWriter->Update();


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

	return EXIT_SUCCESS;
}
