#include <AntTweakBar.h>
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <stdio.h>


#include <iostream>
// #include <wrap/io_trimesh/import_out.h>
#include <vcg/complex/complex.h>
#include <vcg/complex/append.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/create/ball_pivoting.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/stat.h>
#include <wrap/io_trimesh/export_ply.h>
#include <wrap/io_trimesh/import_ply.h>
#include <vcg/complex/algorithms/update/color.h>
#include <vcg/complex/algorithms/update/selection.h>
#include <vcg/math/matrix44.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/gl/trimesh.h>
#include <wrap/gl/shot.h>

#include "freeimage.h"

    
#include <map>
#ifdef UNIX
#include <unistd.h>
#else
#include <direct.h>
#endif
#include <stdio.h>
#include <fstream>



#define SCALE 0.5
#define CREATE_SINGLE_MESH 1
#define CREATE_MLP 0


#include <wrap/gui/trackball.h>
#include <wrap/gl/shot.h>
#include <wrap/gl/camera.h>

#define CACERES_DATASET_BLUE
const bool tessellate = false;

int bin_step = 10;
int img_step = 10;

using namespace std;

class MyVertex; class MyEdge; class MyFace;
struct MyUsedTypes : public vcg::UsedTypes<vcg::Use<MyVertex>   ::AsVertexType,
                                           vcg::Use<MyEdge>     ::AsEdgeType,
                                           vcg::Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public vcg::Vertex< MyUsedTypes, vcg::vertex::BitFlags,vcg::vertex::VFAdj,vcg::vertex::Coord3f, vcg::vertex::Normal3f ,vcg::vertex::Color4b,vcg::vertex::Qualityf>{};
class MyFace    : public vcg::Face<   MyUsedTypes,vcg::face::Mark,vcg::face::VertexRef,vcg::face::BitFlags,vcg::face::Normal3f,vcg::face::FFAdj> {};
class MyEdge    : public vcg::Edge<   MyUsedTypes> {};

class MyMesh    : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> , std::vector<MyEdge>  > {};


vcg::Point3f velodyne_rpy;
vcg::Point3f left_rpy;
vcg::Matrix44d lidar2imu,camera2imu;


std::vector<std::pair<double,vcg::Shotf> > cameras;

// vector of camera extrinsincs
std::vector<vcg::Shotf> shots;

// vector of camera times ( i-th element in image_stamps is the time of i-th camera in shots)
std::vector<double> image_stamps;

// vector of image names
std::vector<std::string> image_names;

// vector of point clouds  names
std::vector<std::string> converted_bins;

struct CameraIntrinsics{
    float focal;
    int c[2];
    int vp[2];
};

CameraIntrinsics camera;

float scale = 1.f;

void import_from_bundle(char*filename)
    {
        MyMesh m;
        MyMesh allCams;
//        vcg::tri::io::ImporterOUT<MyMesh>::Open(m,shots,filename);
        vcg::tri::io::ExporterPLY<MyMesh>::Save(m,"reco.ply",vcg::tri::io::Mask::IOM_ALL);
    }


void readImagesTime(std::string  folder){
    char _[100];
    int s,ms;
   // _chdir(folder.c_str());
    FILE * fi;
    int i=0;
    do{
        fi = fopen((folder+"\\left_"+std::to_string(i)+"_timestamp.txt").c_str(),"r");
        if(fi!=0)
            {
             fgets(_,100,fi);
             fscanf(fi,"%d,%d",&s,&ms);
             fclose(fi);
             image_stamps.push_back(s+0.000000001*double(ms));
             image_names.push_back("left_"+std::to_string(i)+".jpg");
            }
        i += img_step; // this is because we process only image iwth %img_step index
    } while(fi!=0);
}

std::vector< std::pair<double, vcg::Quaterniond>> imu_rotations;
std::vector< std::pair<double, vcg::Point3d>> lidar_translation;
std::vector< std::pair<double, vcg::Quaterniond>> lidar_rotation;




vcg::Quaterniond rotationFromImu(double ts){
    int i = 0;
    do{i++;}
    while (imu_rotations[i].first < ts );
    double alpha = (ts - imu_rotations[i-1].first)/(imu_rotations[i].first - imu_rotations[i-1].first);
    alpha = std::min(1.0,std::max(0.0,alpha));
    return vcg::Interpolate<double>(imu_rotations[i-1].second,imu_rotations[i].second,alpha);
}


vcg::Point3d positionFromBin(double ts){
    int i = 0;
    do{i++;}
    while ( i<lidar_translation.size() && lidar_translation[i].first < ts );
    double alpha = (ts - lidar_translation[i-1].first )/(lidar_translation[i].first  - lidar_translation[i-1].first);
    return lidar_translation[i-1].second* (1.0-alpha) + lidar_translation[i].second * alpha;
}

vcg::Quaterniond rotationFromBin(double ts, int &l){
    int i = 0;
    do{i++;}
    while ( i+1<lidar_rotation.size() && lidar_rotation[i].first < ts );
    l=i-1;
    double alpha = (ts - lidar_rotation[i-1].first )/(lidar_rotation[i].first  - lidar_rotation[i-1].first);
    return vcg::Interpolate<double>(lidar_rotation[i-1].second,lidar_rotation[i].second,alpha);
}

//std::vector< std::pair<double, vcg::Quaterniond>> imu_rotations;

vcg::Shotf cameraFromSfM(double ts){
    int i = 0;
    do{i++;}
    while ( i<image_stamps.size() && image_stamps[i] < ts );

    printf("ts: %d of %d \n",i,image_stamps.size());

    double alpha = (ts - image_stamps[i-1] )/(image_stamps[i]  - image_stamps[i-1]);
    alpha = std::min(1.0,std::max(0.0,alpha));

    vcg::Shotf res;
    res.SetViewPoint( shots[i-1].GetViewPoint()*(1.f-alpha)+shots[i].GetViewPoint()*alpha);
    vcg::Quaternionf q1,q2,q_int;
    q1.FromMatrix(shots[i-1].Extrinsics.Rot());
    q2.FromMatrix(shots[i].Extrinsics.Rot());
    q_int = vcg::Interpolate<float>(q1,q2,alpha);
    vcg::Matrix44f m_int;
    q_int.ToMatrix(m_int);
    res.Extrinsics.SetRot(m_int);
    return res;
}

MyMesh allLerpCams;

void RemoveLongFaces(MyMesh & mm,float perc){

    vcg::Distribution<float> d;

    for(MyMesh::FaceIterator fi = mm.face.begin(); fi != mm.face.end(); ++fi)
        for(unsigned int i = 0;i < 3; ++i)
            d.Add(vcg::Distance((*fi).P(i),(*fi).P((i+1)%3)));

    vcg::tri::UpdateSelection<MyMesh>::FaceOutOfRangeEdge(mm,0,d.Percentile(perc));
     for(MyMesh::FaceIterator fi = mm.face.begin(); fi != mm.face.end(); ++fi)
         if(!(*fi).IsD() && (*fi).IsS())
             vcg::tri::Allocator<MyMesh>::DeleteFace(mm,*fi);

}


template <class MST,class VST>
vcg::Point3d matmul(vcg::Matrix44<MST> mat, vcg::Point3<VST> p){
    vcg::Point3d pf;
    vcg::Matrix44d matf;
    matf.Import(mat);
    pf.Import(p);
    vcg::Point4d p4 = matf*vcg::Point4d(pf[0],pf[1],pf[2],1.0);
    return *(vcg::Point3d*)&p4;
}
vcg::Point2f vcgPoint2f(vcg::Point3f p){return vcg::Point2f(p[0],p[1]);}

int readBinMesh_no_T(std::string meshname){
    printf("processing %s\n",meshname.c_str());

    MyMesh m;

    vcg::Point3d tr;
    double qr[4];

    FILE *f = fopen(meshname.c_str(),"rb");
    if(f==0)
        return (0);
    int n_points;


    int t[2];
    fread(t,sizeof(int),2,f);


    fread(&tr[0],sizeof(double),3,f);
    fread(&qr[0],sizeof(double),4,f);
    fread(&n_points,sizeof(int),1,f);


    MyMesh::VertexIterator vi = vcg::tri::Allocator<MyMesh>::AddVertices(m,n_points);
    for(int i = 0; i < n_points;++i){
        fread(&(*vi).P(),sizeof(float),3,f);
        fread(&(*vi).Q() ,sizeof(float),1,f);
        ++vi;
    }
    vcg::tri::io::ExporterPLY<MyMesh>::Save(m,(meshname+"_noT.ply").c_str(),vcg::tri::io::Mask::IOM_ALL);
return 1;
}


MyMesh total;
int readBinMeshNoSfm(std::string meshname, bool tessellate = true){
 
    char buffer[65536];
    _getcwd(buffer, 1000);
    printf("dir: %s\n", buffer);
   printf("processing %s\n",meshname.c_str());
    MyMesh m;
    int t[2]; // timestamp seconds, nanoseconds
    vcg::Point3d tr;
    vcg::Quaterniond q;
    double qr[4];
    int n_points;


    FILE *f = fopen(meshname.c_str(),"rb");
    if(f==0)
        return (0);



    fread(t,sizeof(int),2,f);
    double ts = t[0]+0.000000001*double(t[1]);

  //  vcg::Quaterniond rotImu = rotationFromImu(ts);

    fread(&tr[0],sizeof(double),3,f);

    fread(&qr[0],sizeof(double),4,f);
    q[0] = qr[3];
    *(vcg::Point3d*)&q[1] = *(vcg::Point3d*)&qr[0];
   
    fread(&n_points,sizeof(int),1,f);

    // ptc writing error. No rotation is not a null quaternion
    bool use_q;
    double no = fabs((*(vcg::Point3d*)&q[1]).Norm()-1.0);
    use_q = no < 0.1;

    MyMesh mor;
    MyMesh::VertexIterator vi = vcg::tri::Allocator<MyMesh>::AddVertices(m,n_points);
    std::vector<vcg::Point3f> pos;
    pos.resize(n_points);

    for(int i = 0; i < n_points;++i){
        float ro,theta,phi;
        fread(&pos[i],sizeof(float),3,f);
        fread(&(*vi).Q() ,sizeof(float),1,f);

        if(tessellate){
            vcg::Point3f p2 = pos[i];
            vcg::Point3f p1;
            p1[0] = p2[0];
            p1[1] = p2[2];
            p1[2] = -p2[1];
            p1.ToPolarRad(ro,theta,phi);
            (*vi).P() = vcg::Point3f(theta,phi,0);
        }else
        (*vi).P() = pos[i];
        ++vi;
        }

    if(tessellate){
        vcg::tri::BallPivoting<MyMesh> bpivot(m,0,0);
        bpivot.BuildMesh();

        for(MyMesh::FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
            if(   ((vcgPoint2f(fi->cP(1))-vcgPoint2f(fi->cP(0))) ^
                  (vcgPoint2f(fi->cP(2))-vcgPoint2f(fi->cP(0))) ) < 0.f)
                    std::swap(fi->V(0),fi->V(1));
        vi =  m.vert.begin();
        for(unsigned int i=0 ; i<pos.size();++i){
           (*vi).P()  = pos[i];
            ++vi;
        }

        vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);
        vcg::tri::Clean<MyMesh>::RemoveDuplicateFace(m);
        vcg::tri::Allocator<MyMesh>::CompactFaceVector(m);
        vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);
        RemoveLongFaces(m,0.35);
        vcg::tri::Allocator<MyMesh>::CompactFaceVector(m);
        vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);
        vcg::tri::Clean<MyMesh>::RemoveSmallConnectedComponentsSize(m,100);
        }
//    vcg::tri::io::ExporterPLY<MyMesh>::Save(m,(meshname+"_T.ply").c_str(),vcg::tri::io::Mask::IOM_ALL);

    for(MyMesh::VertexIterator vi = m.vert.begin(); vi != m.vert.end();++vi){

        vcg::Point3d pd;

        pd.Import(vi->P()*scale);


        if(use_q)
            pd = q.Rotate(pd);

        pd+=tr;

        vi->P().Import(pd);
        }



    // tessellate
    if(tessellate)
        vcg::tri::io::ExporterPLY<MyMesh>::Save(m,(meshname+".T.ply").c_str(),vcg::tri::io::Mask::IOM_ALL);
    else
        vcg::tri::io::ExporterPLY<MyMesh>::Save(m,(meshname+".ply").c_str(),vcg::tri::io::Mask::IOM_ALL);

    fclose(f);

    if(CREATE_SINGLE_MESH)
       vcg::tri::Append<MyMesh, MyMesh>::Mesh(total, m);

    return 1;
}


void readImuTS(std::string filename){
    FILE*f = fopen(filename.c_str(),"rb");
    if(f==0)
        exit(0);
    char _[100];
    fgets(_,100,f);
//    vcg::Matrix44d M;
//    int i = 0;
    while(!feof(f)){
        int s,n;
        float x,y,z,w;
        fscanf(f,"%d, %d, %f, %f, %f, %f",&s,&n,&x,&y,&z,&w);
        imu_rotations.push_back(std::pair<double, vcg::Quaterniond>());
        imu_rotations.back().first=s+0.000000001*double(n);
        imu_rotations.back().second.X() = w;//
        imu_rotations.back().second.Y() = x;//
        imu_rotations.back().second.Z() = y;//
        imu_rotations.back().second.W() = z;//
        imu_rotations.back().second.Normalize();
    }
    fclose(f);
}

int readBinTranslation(std::string meshname){
    int t[2]; // timestamp seconds, nanoseconds
    vcg::Quaterniond q;
    double qr[4];

    vcg::Point3d tr;
    FILE *f = fopen(meshname.c_str(),"rb");
    if(f==0)
        return  0 ;

    fread(t,sizeof(int),2,f);
    double ts = t[0]+0.000000001*double(t[1]);
    fread(&tr[0],sizeof(double),3,f);
    lidar_translation.push_back(std::make_pair(ts,tr));

    fread(&qr[0],sizeof(double),4,f);
    q[0] = qr[3];
    *(vcg::Point3d*)&q[1] = *(vcg::Point3d*)&qr[0];
    lidar_rotation.push_back(std::make_pair(ts,q));

    fclose(f);
    return 1;
}

int readBinMesh(std::string meshname){
    printf("processing %s\n",meshname.c_str());
    MyMesh m;
    int t[2]; // timestamp seconds, nanoseconds
    vcg::Point3d tr;
    vcg::Quaterniond q;
    double qr[4];
    int n_points;

    FILE *f = fopen(meshname.c_str(),"rb");
    if(f==0)
        return (0);

    fread(t,sizeof(int),2,f);
    double ts = t[0]+0.000000001*double(t[1]);

    vcg::Quaterniond imuQ = rotationFromImu(ts);

    fread(&tr[0],sizeof(double),3,f);
    fread(&qr[0],sizeof(double),4,f);
    *(vcg::Point3d*)&q[1] = *(vcg::Point3d*)&qr[0];
    q[0] = qr[3];
    fread(&n_points,sizeof(int),1,f);

    vcg::Shotf shot = cameraFromSfM(ts);
    tr.Import(shot.GetViewPoint());

    // check for writing error. No rotation is not a null quaternion
    bool use_q;
    use_q = fabs((*(vcg::Point3d*)&q[1]).Norm()-1.0) < 0.01;

    MyMesh mor;
    MyMesh::VertexIterator vi = vcg::tri::Allocator<MyMesh>::AddVertices(m,n_points);
    std::vector<vcg::Point3f> pos;
    pos.resize(n_points);

    for(int i = 0; i < n_points;++i){
        float ro,theta,phi;
        fread(&pos[i],sizeof(float),3,f);
        fread(&(*vi).Q() ,sizeof(float),1,f);

        vcg::Point3f p2 = pos[i];
        vcg::Point3f p1;
        p1[0] = p2[0];
        p1[1] = p2[2];
        p1[2] = -p2[1];

        p1.ToPolarRad(ro,theta,phi);

        (*vi).P() = vcg::Point3f(theta,phi,0);
        ++vi;
        }
    vcg::tri::BallPivoting<MyMesh> bpivot(m,0,0);
    bpivot.BuildMesh();

    for(MyMesh::FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
        if(   ((vcgPoint2f(fi->cP(1))-vcgPoint2f(fi->cP(0))) ^
              (vcgPoint2f(fi->cP(2))-vcgPoint2f(fi->cP(0))) ) < 0.f)
                std::swap(fi->V(0),fi->V(1));

    vi =  m.vert.begin();
    for(unsigned int i=0 ; i<pos.size();++i){
       (*vi).P()  = pos[i];
        ++vi;
    }

    vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);
    vcg::tri::Clean<MyMesh>::RemoveDuplicateFace(m);
    vcg::tri::Allocator<MyMesh>::CompactFaceVector(m);
    vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);
    RemoveLongFaces(m,0.35);
    vcg::tri::Allocator<MyMesh>::CompactFaceVector(m);
    vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);
    vcg::tri::Clean<MyMesh>::RemoveSmallConnectedComponentsSize(m,100);

//    vcg::tri::io::ExporterPLY<MyMesh>::Save(m,(meshname+"_T.ply").c_str(),vcg::tri::io::Mask::IOM_ALL);

    for(MyMesh::VertexIterator vi = m.vert.begin(); vi != m.vert.end();++vi){

        vcg::Point3d pd;

        pd.Import(vi->P()*scale);
        if(use_q)
            pd = q.Rotate(pd);

        // apply lidar2imuransform to bring the point cloud to the velodyne frame
       pd = matmul(lidar2imu,pd);

        //pd = imuQ.Rotate(pd);

        // apply inverse of camera2imu to bring the point cloud to the camera frame (mounted on the velodyne)
       // pd = matmul(camera2imu.transpose(),pd);
       pd = matmul(vcg::Inverse(camera2imu),pd);

        // conversion ROS to VCG
        vcg::Matrix44f mat;
        pd = matmul(mat.SetRotateDeg(180,vcg::Point3f(1,0,0)),pd);

        // project the point in the global frame
        vcg::Point3f pdf;
        pd = matmul(shot.Extrinsics.Rot().transpose(),pd)+tr;

        vi->P().Import(pd);
        }

    // tessellate

    vcg::tri::io::ExporterPLY<MyMesh>::Save(m,(meshname+".ply").c_str(),vcg::tri::io::Mask::IOM_ALL);

    fclose(f);

    return 1;
}

void align_bins(std::string folder){
    _chdir(folder.c_str());
    int i=0;
    while(readBinMesh((folder+"\\"+std::to_string(i)+".bin"))){i+=bin_step;}
}

vcg::Matrix44f rpy2mat(vcg::Point3f rpy){
    vcg::Matrix44f res,rot;
    res.SetRotateRad(rpy[0],vcg::Point3f(1,0,0));
    rot.SetRotateRad(rpy[1],vcg::Point3f(0,1,0));    res = rot * res;
    rot.SetRotateRad(rpy[2],vcg::Point3f(0,0,1));
    res = rot * res;
    return res;
}



void readLidarTranslations(char * folder){
    int i=0;
    while(readBinTranslation((std::string(folder)+"\\"+std::to_string(i)+".bin"))){
        i++;
    }
}

void bins2ply(char * folder){
    int i=0;
    while( readBinMeshNoSfm((std::string(folder)+'\\'+std::to_string(i)+".bin"),tessellate)){
        i+=bin_step;
        converted_bins.push_back(std::to_string(i)+".bin");
    }
    if (CREATE_SINGLE_MESH)
        vcg::tri::io::ExporterPLY<MyMesh>::Save(total, "total.ply");
}


#ifdef CACERES_DATASET_BLUE
double T_cam_imu[16]={
  -0.026367054322066402, -0.9996271606687501, 0.00709352519636397, 0.023448320704025388,
  -0.030838652538258182, -0.006279228147363675, -0.9995046517167889, 0.02906892970104672,
  0.999176538933938, -0.026572748205776653, -0.03066158992600701, -0.14528314791638972,
  0.0, 0.0, 0.0, 1.0
};

double T_scan_imu[16]={
    0.9535461825, 0, -0.3012468719, -0.08483511244,
    0,          -1.0, 0,                0.01,
    -0.3012468719,0, -0.9535461825,  -0.1357497834,
    0,    0,             0,            1.0,
    };
#endif


vcg::Point3d skew(vcg::Matrix44d m){
    return vcg::Point3d(m.GetColumn3(0)*m.GetColumn3(1),
    m.GetColumn3(0)*m.GetColumn3(2),
    m.GetColumn3(1)*m.GetColumn3(2));
}
 

vcg::Matrix44d lidarAtTime(double t) {
    vcg::Matrix44d RR;

    vcg::Point3d p = positionFromBin(t);
    int il;
    vcg::Quaterniond qrotL = rotationFromBin(t, il);
    vcg::Matrix44d lidarFrame;
    qrotL.ToMatrix(lidarFrame);

    lidarFrame.SetColumn(3, p);

    return lidarFrame;
}


vcg::Shotd  cameraAtTime(double t) {
    vcg::Matrix44d RR;
    vcg::Matrix44d cameraFrame;
    vcg::Matrix44d imu2camera = vcg::Inverse(camera2imu);
    vcg::Matrix44d lidarFrame;
    vcg::Point3d p;
    vcg::Shotd shot;

    int il;
    p = positionFromBin(t);
    vcg::Quaterniond qrotL = rotationFromBin(t, il);
    qrotL.ToMatrix(lidarFrame);
    lidarFrame.SetColumn(3, p);
    cameraFrame = lidarFrame * lidar2imu  * imu2camera;
   
    shot.SetViewPoint(cameraFrame.GetColumn3(3));
    cameraFrame.SetColumn(3, vcg::Point4d(0, 0, 0, 1));
    cameraFrame = cameraFrame * RR.SetRotateDeg(180, vcg::Point3d(1, 0, 0));
    shot.Extrinsics.SetRot(cameraFrame);
    shot.Intrinsics.SetFrustum(-camera.vp[0] / 2, camera.vp[0] / 2, -camera.vp[1] / 2, camera.vp[1] / 2, camera.focal, vcg::Point2i(camera.vp[0], camera.vp[1]));

    return shot;
}
 
vcg::Matrix44d  imuAtTime(double t) {
    vcg::Matrix44d imuFrame;
    vcg::Matrix44d RR;
    vcg::Matrix44d cameraFrame;
    vcg::Matrix44d imu2camera = vcg::Inverse(camera2imu);
    vcg::Matrix44d lidarFrame;
    vcg::Point3d p;
    vcg::Shotd shot;

    int il;
    p = positionFromBin(t);
    vcg::Quaterniond qrotL = rotationFromBin(t, il);
    qrotL.ToMatrix(lidarFrame);
    lidarFrame.SetColumn(3, p);
    imuFrame = lidarFrame * lidar2imu ;

    return imuFrame;
}

void create_mlp(char * folder, bool tessellate ){
   std::string fname = std::string(folder)+"/../meshlab.mlp";

    vcg::Matrix44d RR;
    vcg::Matrix44d rot;
    vcg::Matrix44d imuFrame;
    vcg::Matrix44d cameraFrame;
    vcg::Matrix44d imu2camera;

    std::string mname = (tessellate)?"_T.ply":".ply";
    vcg::Point3d p;
    std::ofstream of;
    of.open(fname);
    of <<"<!DOCTYPE MeshLabDocument>\n <MeshLabProject> "<< std::endl;

    of <<"<MeshGroup>" << std::endl;
    for(int i=0; i < converted_bins.size()-1; i++)
    of <<"<MLMesh label=\""<< converted_bins[i] <<"\" filename=\""<< "LIDAR/"<<converted_bins[i]<<mname <<"\">\n<MLMatrix44> 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 </MLMatrix44> </MLMesh>" << std::endl;
    of <<"</MeshGroup>\n";

    of << "<RasterGroup>" << std::endl;
    for(int i=0; i <  image_names.size();i++){

        vcg::Shotd shot = cameraAtTime(image_stamps[i]);
        rot = shot.Extrinsics.Rot().transpose();
        p = shot.GetViewPoint();
   
    of << "<MLRaster label=\""  << image_names[i].c_str() << "\">\n"
       << "<VCGCamera FocalMm=\""<< camera.focal << "\" RotationMatrix=\""<<
          rot[0][0] << " " << rot[0][1] <<" "<< rot[0][2] << " " << 0<< " "<<
          rot[1][0] << " " << rot[1][1] <<" "<< rot[1][2] << " " << 0<< " "<<
          rot[2][0] << " " << rot[2][1] <<" "<< rot[2][2] << " " << 0<< " "<<
          rot[3][0] << " " << rot[3][1] <<" "<< rot[3][2] << " " << 1<<
          "\" ViewportPx=\""<< camera.vp[0] << " " << camera.vp[1] <<"\" LensDistortion=\"0 0\"  BinaryData=\"0\"  CameraType=\"0\" TranslationVector=\""<<  -p[0] << " " << -p[1] << " "<<-p[2] <<" 1\" PixelSizeMm=\""<< 1.f   <<" " << 1.f <<"\""<<std::endl;
    of<< "CenterPx=\""<< camera.c[0]<<" " << camera.c[1]<<"\"/>\n"<<
          "<Plane semantic=\"1\" fileName=\""<< "camera/"<<image_names[i].c_str() <<"\"/>\n</MLRaster>"<<std::endl;
    }
    of << " </RasterGroup></MeshLabProject>" << std::endl;
    of.close();
}

void stats_lidar(char * folder){
    vcg::Quaterniond q1,q2,q;
    std::string fname = std::string(folder)+"/../rots.txt";
    FILE*f=fopen(fname.c_str(),"w");
    double phi,t;
    vcg::Point3d _,p1,p2;
    for(int i= 0; i < lidar_rotation.size()-1; ++i){
        t = lidar_rotation[i+1].first-lidar_rotation[i].first;

        q1 = lidar_rotation[i].second;
        q2  = lidar_rotation[i+1].second;
         q1.Invert();
        q =q1*q2;
        q.ToAxis(phi,_);
        p1 = lidar_translation[i].second;
        p2 = lidar_translation[i+1].second;

        fprintf(f,"%d %f %f\n",i,phi/t,(p2-p1).Norm()/t);
    }
    fclose(f);
}


void print_usage(){
    std::cout << "dron2home PATH_TO_imu.txt PATH_TO_IMAGES_FOLDER PATH_TO_bin_FOLDER bin_step{default=10} img_step{default=10} "<< std::endl;
    std::cout << "It takes the data collected from the drone and it creates a meshlab project\n\
 bin_step and img_step can be used to choose a subset of the point clouds and images\n\
 You can set them to 1 and the whole data will be processed but then you want be able to open\n\
 the mlp file with meshlab" << std::endl;
}
std::string images_folder;

int mainGUI(int argc, char** argv);
int main(int argc,char ** argv)
{
    
//    readBinMeshNoSfm((std::string(argv[3])+"/"+std::to_string(1640)+".bin"),true);
//    return 0;

    // dataset july
    // velodyne_rpy = vcg::Point3f(3.151592653589793,0.1323284641020683,0.);
    // left_rpy = vcg::Point3f(-1.5708, 0, -1.5708);

    // lidar2imu = rpy2mat(velodyne_rpy);
    // camera2imu = rpy2mat(left_rpy);


    print_usage();

#ifdef CACERES_DATASET_BLUE
    camera2imu = vcg::Matrix44d(T_cam_imu);

    lidar2imu = vcg::Matrix44d(T_scan_imu);



    camera.focal = 1400.71;
    camera.c[0]  = 1017.87000;
    camera.c[1]  = 553.07000;
    camera.vp[0] = 1920;
    camera.vp[1] = 1080;

#endif

    scale = 1.0;
    if(argc>4){
        bin_step = atoi(argv[4]);
        img_step = atoi(argv[5]);
    }

    readLidarTranslations(argv[3]);
    // read imu values
    readImuTS(std::string(argv[1]));

    images_folder = std::string(argv[2]);
    // read image time stamps
    readImagesTime(std::string(argv[2]));

   
    //bins2ply(argv[3]);

    if (CREATE_MLP) {
        // convert bin to plys
        bins2ply(argv[3]);
        // create MLP
        create_mlp(argv[3], tessellate);
        stats_lidar(argv[3]);
    }

    mainGUI(  argc,   argv);
    glutInit(&argc, argv);
}


void Reshape(int,int){}
void terminate(){}

/// we choosed a subset of the avaible drawing modes
enum DrawMode{SMOOTH=0,PERPOINTS,WIRE,FLATWIRE,HIDDEN,FLAT};

/// the current drawmode
DrawMode drawmode;

void TW_CALL CopyCDStringToClient(char **destPtr, const char *src)
{
    size_t srcLen = (src!=NULL) ? strlen(src) : 0;
    size_t destLen = (*destPtr!=NULL) ? strlen(*destPtr) : 0;

    // Alloc or realloc dest memory block if needed
    if( *destPtr==NULL )
        *destPtr = (char *)malloc(srcLen+1);
    else if( srcLen>destLen )
        *destPtr = (char *)realloc(*destPtr, srcLen+1);

    // Copy src
    if( srcLen>0 )
        strncpy(*destPtr, src, srcLen);
    (*destPtr)[srcLen] = '\0'; // null-terminated string
}




// GRAPHICS DEBUG
int width, height;

MyMesh m;

using namespace vcg;

vcg::Trackball track;
vcg::GlTrimesh<MyMesh> glWrap;

GLuint vBuffer;

void initForGUI() {
    vcg::tri::io::ImporterPLY<MyMesh>::Open(m, "all_geometry.ply");
    glWrap.m = &m;
    glWrap.Update();
    vcg::tri::UpdateBounding<MyMesh>::Box(m);

    GLuint id_tex;
    glGenTextures(1, &id_tex);
    glBindTexture(GL_TEXTURE_2D, id_tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1920,1080, 0, GL_BGR, GL_UNSIGNED_BYTE, NULL);
}

static vcg::Trackball::Button GLUT2VCG(int glut_button, int)
{
    int vcgbt = vcg::Trackball::BUTTON_NONE;

    switch (glut_button) {
    case GLUT_LEFT_BUTTON: vcgbt |= vcg::Trackball::BUTTON_LEFT;	break;
    case GLUT_MIDDLE_BUTTON: vcgbt |= vcg::Trackball::BUTTON_RIGHT;	break;
    case GLUT_RIGHT_BUTTON: vcgbt |= vcg::Trackball::BUTTON_MIDDLE;	break;
    }

    int modifiers = glutGetModifiers();

    if (modifiers & GLUT_ACTIVE_SHIFT)	vcgbt |= vcg::Trackball::KEY_SHIFT;
    if (modifiers & GLUT_ACTIVE_CTRL)	vcgbt |= vcg::Trackball::KEY_CTRL;
    if (modifiers & GLUT_ACTIVE_ALT)	vcgbt |= vcg::Trackball::KEY_ALT;

    return vcg::Trackball::Button(vcgbt);
}

void   keyReleaseEvent(unsigned char k, int x, int y)
{
    int modifiers = glutGetModifiers();
    if (modifiers & GLUT_ACTIVE_CTRL)
        track.ButtonUp(Trackball::Button::KEY_CTRL);
    if (modifiers & GLUT_ACTIVE_SHIFT)
        track.ButtonUp(Trackball::Button::KEY_SHIFT);
    if (modifiers & GLUT_ACTIVE_ALT)
        track.ButtonUp(Trackball::Button::KEY_ALT);
}

void   keyPressEvent(unsigned char k, int x, int  y)
{

    int modifiers = glutGetModifiers();
    if (modifiers & GLUT_ACTIVE_CTRL)
        track.ButtonDown(Trackball::Button::KEY_CTRL);
    if (modifiers & GLUT_ACTIVE_SHIFT)
        track.ButtonDown(Trackball::Button::KEY_SHIFT);
    if (modifiers & GLUT_ACTIVE_ALT)
        track.ButtonDown(Trackball::Button::KEY_ALT);

    TwEventKeyboardGLUT(k, x, y);
}

void mousePressEvent(int bt, int state, int x, int y) {
    if (TwEventMouseButtonGLUT(bt, state, x, y))
        return;


        if (state == GLUT_DOWN)
            track.MouseDown(x, height - y, GLUT2VCG(bt, state));
        else
            track.MouseUp(x, height - y, GLUT2VCG(bt, state));


};

void mouseMoveEvent(int x, int y)
{
        if (!TwEventMouseMotionGLUT(x, y))
            track.MouseMove(x, height - y);
}



void wheelEvent(int wheel, int direction, int x, int y) {
    track.MouseWheel(wheel * direction);
}


void* currentImage;
double currentTime;
int n_image=-1;
vcg::Shotd  currentShot;
vcg::Matrix44d currentLidar;
vcg::Matrix44d currentImu;

double near_plane=1.0, far_plane=100.0;

void draw_axes(vcg::Point3f col) {

    glLineWidth(3.0);
    glBegin(GL_LINES);
    glColor3f(col[0], 0, 0);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(1.0, 0.0, 0.0);


    glColor3f(0, col[1], 0);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 1.0, 0.0);

    glColor3f(0, 0, col[2]);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 0.0, 1.0);
    glEnd();

}
void draw_shot(vcg::Shotd s) {
    vcg::Point3d p = s.GetViewPoint();
    glPushMatrix();
    glTranslated(p[0], p[1], p[2]);
    vcg::Matrix44d R = s.Extrinsics.Rot().transpose();
    glMultMatrixd(&R[0][0]);
    draw_axes(vcg::Point3f(1, 1, 1));
    glPopMatrix();

}
void draw_frame(vcg::Matrix44d frame) {

    glPushMatrix();
    vcg::Matrix44d R = frame.transpose();
    glMultMatrixd(&R[0][0]);
    draw_axes(vcg::Point3f(0.5, 0.5, 0.5));
    glPopMatrix();
}

void print_matrix(const char* n, double * m) {
    printf("%s\n %f %f %f %f \n %f %f %f %f \n %f %f %f %f \n %f %f %f %f\n", n, m[0], m[1], m[2], m[3],
        m[4], m[5], m[6], m[7],
        m[8], m[9], m[10], m[11],
        m[12], m[13], m[14], m[15]);
}


void Display() {

    glClearColor(0.0, 0.32, 0.48, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDisable(GL_LIGHTING);

    if(1) {
        glDepthRange(0.5, 1.0);
        glViewport(width/2 , 0, width/2,height/2 );
        //MONITOR VIEW
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glFrustum(-1920 / 2 * 0.001, 1920 / 2 * 0.001, -1080 / 2 * 0.001, 1080 / 2 * 0.001, 1.4, far_plane);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        gluLookAt(0, 0, 5, 0, 0, 0, 0, 1, 0);

        track.center = vcg::Point3f(0, 0, 0);
        track.radius = 1;
        track.GetView();
        track.Apply();
        glPushMatrix();
        float d = 1.0f / m.bbox.Diag();
        vcg::glScale(d);
        glTranslate(-glWrap.m->bbox.Center());
        if (0) {
            glPushMatrix();
            glScaled(10., 10., 10.);
            draw_axes(vcg::Point3d(1, 1, 1));
            glPopMatrix();
        }

        glWrap.Draw<vcg::GLW::DMPoints, vcg::GLW::CMPerVert, vcg::GLW::TMNone>();
        if (currentShot.Intrinsics.FocalMm > 0.0) {
            draw_shot(currentShot);
            draw_frame(currentLidar);
            draw_frame(currentImu);
        }
        glPopMatrix();
        track.DrawPostApply();
    }

    if(currentShot.Intrinsics.FocalMm>0.0)
    {
        glDepthRange(0.0, 0.5);

        // from camera view
        glDisable(GL_LIGHTING);
        glColor3f(1.0, 1.0, 1.0);
        glViewport(0, 0, 1920 / 2* SCALE, 1080 / 2* SCALE);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        assert(glGetError() == 0);
          
        glDepthMask(GL_FALSE);
        glEnable(GL_TEXTURE_2D);
//        glEnable(GL_BLEND);
 //       glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        assert(glGetError() == 0);

        glBegin(GL_QUADS);
        glTexCoord2f(0.0, 0.0);
        glVertex3f(-1, -1, 0);
        glTexCoord2f(1.0, 0.0);
        glVertex3f(1, -1, 0);
        glTexCoord2f(1.0, 1.0);
        glVertex3f(1, 1, 0);
        glTexCoord2f(0.0, 1.0);
        glVertex3f(-1, 1, 0);
        glEnd();

 //       glDisable(GL_BLEND);
        glDisable(GL_TEXTURE_2D);
         
        if (far_plane < near_plane)
            far_plane = near_plane + 1.0;
        glDepthMask(GL_TRUE);

 
        vcg::Shotd _s = currentShot;
        _s.Extrinsics.SetRot(_s.Extrinsics.Rot().transpose());

        GlShot<vcg::Shotd>::SetView(_s, 1.4, far_plane);

        glWrap.Draw<vcg::GLW::DMPoints, vcg::GLW::CMPerVert, vcg::GLW::TMNone>();
        GlShot<vcg::Shotd>::UnsetView();


    }


    TwDraw();

    // Present frame buffer
    glutSwapBuffers();

    // Recall Display at next frame
    glutPostRedisplay();

}

int time_shift = 0;
void TW_CALL setShift(const void* value, void* clientData) {
    time_shift = *(int*)value;
    currentShot = cameraAtTime(::image_stamps[n_image]+ double(time_shift)*0.01);
    currentLidar = lidarAtTime(::image_stamps[n_image] + double(time_shift) * 0.01);
    currentImu = imuAtTime(::image_stamps[n_image] + double(time_shift) * 0.01);
}

void TW_CALL getShift(void* value, void* clientData) 
{
    *(int*)value = time_shift;  // for instance
}

void TW_CALL setImage(const void* value, void* clientData)
{

    char buffer[65536];
    _getcwd(buffer, 1000);
    printf("dir: %s\n", buffer);


    n_image = *(const int*)value;  // for instance
    currentTime = ::image_stamps[n_image];
    currentShot = cameraAtTime(::image_stamps[n_image]);
 

    currentLidar = lidarAtTime(::image_stamps[n_image]);
    currentImu = imuAtTime(::image_stamps[n_image]);

    std::string imPath = images_folder + "\\"+image_names[n_image];
    FIBITMAP* bitmap = FreeImage_Load(FIF_JPEG, imPath.c_str(), JPEG_DEFAULT);
    
    BYTE * image_data = FreeImage_GetBits(bitmap);
    //get the image width and height
    int image_width = FreeImage_GetWidth(bitmap);
    int image_height = FreeImage_GetHeight(bitmap);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1920, 1080, 0, GL_BGR, GL_UNSIGNED_BYTE, image_data);

    sprintf((char*)currentImage, "%s", ::image_names[n_image].c_str());
  //  time_shift = 0.0;
}
void TW_CALL getImage(void* value, void* clientData)
{
    *(int*)value = n_image;  // for instance
}
int mainGUI(int argc, char** argv){
    FreeImage_Initialise(true);
    



    currentImage = new char[100];
    sprintf((char*)currentImage, "cio");
    ((char*)currentImage)[4] = '\0';

    width  = 2048 * SCALE;
    height = 1152 * SCALE;

  

    TwBar *bar; // Pointer to the tweak bar

    // Initialize AntTweakBar
    // (note that AntTweakBar could also be intialized after GLUT, no matter)
    if( !TwInit(TW_OPENGL, NULL) )
    {
        // A fatal error occured
        fprintf(stderr, "AntTweakBar initialization failed: %s\n", TwGetLastError());
        return 1;
    }

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(width, height);
    glutCreateWindow("AntTweakBar simple example using GLUT");
    glutCreateMenu(NULL);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHTING);

    // Set GLUT callbacks
    glutDisplayFunc(Display);
    glutReshapeFunc(Reshape);
//    atexit(Terminate);  // Called after glutMainLoop ends

        // Set GLUT event callbacks
    // - Directly redirect GLUT mouse button events to AntTweakBar
    glutMouseFunc((GLUTmousebuttonfun)mousePressEvent);
    // - Directly redirect GLUT mouse motion events to AntTweakBar
    glutMotionFunc((GLUTmousemotionfun)mouseMoveEvent);
    // - Directly redirect GLUT mouse "passive" motion events to AntTweakBar (same as MouseMotion)
    glutPassiveMotionFunc((GLUTmousemotionfun)TwEventMouseMotionGLUT);
    // - Directly redirect GLUT key events to AntTweakBar
    glutKeyboardFunc((GLUTkeyboardfun)TwEventKeyboardGLUT);
    // - Directly redirect GLUT special key events to AntTweakBar
    glutSpecialFunc((GLUTspecialfun)TwEventSpecialGLUT);

    glutKeyboardFunc(keyPressEvent);
    glutKeyboardUpFunc(keyReleaseEvent);


   glutMouseWheelFunc(wheelEvent);
    TwWindowSize(width, height);
    bar = TwNewBar("Controls");
    TwDefine("Controls size='240 240'"); // resize bar

    TwCopyCDStringToClientFunc (CopyCDStringToClient);

    TwAddVarRO(bar, "image", TW_TYPE_CSSTRING(100), currentImage, " ");
    TwAddVarRO(bar, "time", TW_TYPE_DOUBLE, &currentTime, "precision=3 ");


    
    TwAddVarCB(bar, "setimage", TW_TYPE_INT32, setImage, getImage, 0, (std::string("min='0' max='")+ std::to_string(image_names.size()-1)+ "' ").c_str());
    TwAddVarCB(bar, "shift", TW_TYPE_INT32, setShift, getShift, 0, "");

    TwAddVarRW(bar,"near",TW_TYPE_DOUBLE,&near_plane," label='near' step='0.1'   ");
    TwAddVarRW(bar,"far",TW_TYPE_DOUBLE,&far_plane," label='far' step='0.1'  ");


  //  TwAddVarRW(bar,"Input",TW_TYPE_CDSTRING,&filename," label='Filepath' group=SetMesh help=` Name of the file to load` ");
  //  TwAddButton(bar,"Load from file",loadMesh,0,	" label='Load Mesh' group=SetMesh help=`load the mesh` ");
  //  TwAddButton(bar,"Use tetrahedron",loadTetrahedron,0,	" label='Make Tetrahedron' group=SetMesh help=`use tetrahedron.` ");
  //  TwAddButton(bar,"Use dodecahedron",loadDodecahedron,0,	" label='Make Dodecahedron' group=SetMesh help=`use dodecahedron.` ");


    // ShapeEV associates Shape enum values with labels that will be displayed instead of enum values
    TwEnumVal drawmodes[6] = { {SMOOTH, "Smooth"}, {PERPOINTS, "Per Points"}, {WIRE, "Wire"}, {FLATWIRE, "FlatWire"},{HIDDEN, "Hidden"},{FLAT, "Flat"}};
    // Create a type for the enum shapeEV
  //  TwType drawMode = TwDefineEnum("DrawMode", drawmodes, 6);
    // add 'g_CurrentShape' to 'bar': this is a variable of type ShapeType. Its key shortcuts are [<] and [>].
  //  TwAddVarRW(bar, "Draw Mode", drawMode, &drawmode, " keyIncr='<' keyDecr='>' help='Change draw mode.' ");

    glewInit();
    initForGUI();

    glutMainLoop();
    return 0;
}
