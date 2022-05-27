#include <iostream>
#include <wrap/io_trimesh/import_out.h>
#include <vcg/complex/complex.h>
#include <vcg/complex/append.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/create/ball_pivoting.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/stat.h>
#include <wrap/io_trimesh/export_ply.h>
#include <vcg/complex/algorithms/update/color.h>
#include <vcg/complex/algorithms/update/selection.h>
#include <vcg/math/matrix44.h>

#include <map>
#include <unistd.h>
#include <stdio.h>
#include <fstream>

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

// vector of image name
std::vector<std::string> image_names;

// vector of point clouds  name
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
        vcg::tri::io::ImporterOUT<MyMesh>::Open(m,shots,filename);
//        vcg::tri::io::ExporterPLY<MyMesh>::Save(m,"bundle.ply",vcg::tri::io::Mask::IOM_VERTCOORD|vcg::tri::io::Mask::IOM_VERTCOLOR);
//        for(int i=0; i < shots.size();++i){
//            m.Clear();
//            vcg::tri::Hexahedron(m);
//            vcg::tri::UpdatePosition<MyMesh>::Scale(m,0.001);
//            vcg::tri::UpdatePosition<MyMesh>::Translate(m,shots[i].GetViewPoint());
//            vcg::tri::UpdateColor<MyMesh>::PerVertexConstant(m,vcg::Color4b::Red);
//            vcg::tri::Append<MyMesh,MyMesh>::Mesh(allCams,m);
//        }
//        vcg::tri::io::ExporterPLY<MyMesh>::Save(allCams,"cams_sfm.ply",vcg::tri::io::Mask::IOM_ALL);
        vcg::tri::io::ExporterPLY<MyMesh>::Save(m,"reco.ply",vcg::tri::io::Mask::IOM_ALL);
    }


void readImagesTime(std::string  folder){
    char _[100];
    int s,ms;
    chdir(folder.c_str());
    FILE * fi;
    int i=0;
    do{
        fi = fopen((folder+"/left_"+std::to_string(i)+"_timestamp.txt").c_str(),"r");
        i += img_step; // this is because we process only image iwth %img_step index
        if(fi!=0)
            {
             fgets(_,100,fi);
             fscanf(fi,"%d,%d",&s,&ms);
             fclose(fi);
             image_stamps.push_back(s+0.000000001*double(ms));
             image_names.push_back("left_"+std::to_string(i)+".jpg");
            }
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
    alpha = std::min(1.0,std::max(0.0,alpha));
    return lidar_translation[i-1].second* (1.0-alpha) + lidar_translation[i].second * alpha;
}

vcg::Quaterniond rotationFromBin(double ts, int &l){
    int i = 0;
    do{i++;}
    while ( i<lidar_rotation.size() && lidar_rotation[i].first < ts );
    l=i-1;
    double alpha = (ts - lidar_rotation[i-1].first )/(lidar_rotation[i].first  - lidar_rotation[i-1].first);
    printf("alpha %f\n",alpha);
    alpha = std::min(1.0,std::max(0.0,alpha));
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



int readBinMeshNoSfm(std::string meshname, bool tessellate = true){
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

    vcg::Quaterniond rotImu = rotationFromImu(ts);

    fread(&tr[0],sizeof(double),3,f);

    fread(&qr[0],sizeof(double),4,f);
    *(vcg::Point3d*)&q[1] = *(vcg::Point3d*)&qr[0];
    q[0] = qr[3];
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


    vi = vcg::tri::Allocator<MyMesh>::AddVertices(m,600);
    for(int i=0; i< 200;++i)
        for(int j = 0; j < 3; ++j){
        vi->P() = vcg::Point3f( 0.0,0.0,0.0);
        vi->P()[j]+= i*0.005;
        vi->C() = vcg::Color4b(0,0,0,255);
        vi->C()[j] = 255;
        ++vi;
     }

    for(MyMesh::VertexIterator vi = m.vert.begin(); vi != m.vert.end();++vi){

        vcg::Point3d pd;

        pd.Import(vi->P()*scale);


        if(use_q)
            pd = q.Rotate(pd);

        pd+=tr;

        vi->P().Import(pd);
        }



//    vcg::Matrix44f x90;
//    x90.SetRotateDeg(90,vcg::Point3f(1,0,0));
//    for(MyMesh::VertexIterator vi = m.vert.begin(); vi != m.vert.end();++vi)
//        (*vi).P() = x90*vi->P();

//        vcg::Matrix44f y180;
//        y180.SetRotateDeg(180,vcg::Point3f(0,1,0));
//        for(MyMesh::VertexIterator vi = m.vert.begin(); vi != m.vert.end();++vi)
//            (*vi).P() = y180*vi->P();


//       // PATCH for MEshlab which goes nuts to up camera vector close to Z
        vcg::Matrix44f x90;
        x90.SetRotateDeg(90,vcg::Point3f(1,0,0));
       for(MyMesh::VertexIterator vi = m.vert.begin(); vi != m.vert.end();++vi)
            (*vi).P() = x90*vi->P();

    // tessellate
    if(tessellate)
        vcg::tri::io::ExporterPLY<MyMesh>::Save(m,(meshname+".T.ply").c_str(),vcg::tri::io::Mask::IOM_ALL);
    else
        vcg::tri::io::ExporterPLY<MyMesh>::Save(m,(meshname+".ply").c_str(),vcg::tri::io::Mask::IOM_ALL);

    fclose(f);
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

//        vcg::Quaterniond q(sqrt(2.0),sqrt(2.0),0.0,0.0);
//        vcg::Quaterniond qn= q * imu_rotations.back().second;

//        vcg::Point3d v,vq,vn;
//        double t,tq,tn;
//        imu_rotations.back().second.ToAxis(t,v);
//        q.ToAxis(tq,vq);
//        qn.ToAxis(tn,vn);


//        printf(" theta,v: %f  (%f, %f, %f), %f \n",t,v[0],v[1],v[2]);
//        printf(" thetaQ,v: %f  (%f, %f, %f), %f \n",tq,vq[0],vq[1],vq[2]);
//        printf(" thetaN,v: %f  (%f, %f, %f), %f \n",tn,vn[0],vn[1],vn[2]);


//        float sinTheta = sin(theta/2.0);
//        vcg::Point3d v(x/sinTheta,y/sinTheta,z/sinTheta);
//        printf(" theta,v: %f  (%f, %f, %f), %f \n",theta,v[0],v[1],v[2],v.Norm());

//        vcg::Quaterniond q1(w,x,y,z);
//        vcg::Quaterniond q2 = q1.Inverse()*imu_rotations.back().second;

//        vcg::Matrix44d m ;
//        q2.ToMatrix(m);

//        vcg::Matrix44d pM;
//        pM.SetZero();
//        pM[0][2] = 1.0;
//        pM[1][0] = -1.0;
//        pM[2][1] = -1.0;
//        pM[3][3] =  1.0;

//        vcg::Quaterniond qM;
//        qM.FromMatrix(pM);

//        vcg::Quaterniond q3 = q1*qM;

//        imu_rotations.back().second = q3;

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
    *(vcg::Point3d*)&q[1] = *(vcg::Point3d*)&qr[0];
    q[0] = qr[3];
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

// saving  lidar positoin as small cubes
//    MyMesh cube;
//    vcg::tri::Hexahedron(cube);
//    vcg::tri::UpdatePosition<MyMesh>::Scale(cube,0.001);
//    vcg::Point3f trf;
//    trf.Import(tr);
//    vcg::tri::UpdatePosition<MyMesh>::Translate(cube,trf);
//    vcg::tri::Append<MyMesh,MyMesh>::Mesh(allLerpCams,cube);
//    vcg::tri::io::ExporterPLY<MyMesh>::Save(cube,(meshname +"_POS_.ply").c_str());

    fclose(f);

    return 1;
}

void align_bins(std::string folder){
    chdir(folder.c_str());
    int i=0;
    while(readBinMesh((folder+"/"+std::to_string(i)+".bin"))){i+=bin_step;}
}

vcg::Matrix44f rpy2mat(vcg::Point3f rpy){
    vcg::Matrix44f res,rot;
    res.SetRotateRad(rpy[0],vcg::Point3f(1,0,0));
    rot.SetRotateRad(rpy[1],vcg::Point3f(0,1,0));    res = rot * res;
    rot.SetRotateRad(rpy[2],vcg::Point3f(0,0,1));
    res = rot * res;
    return res;
}



//int main(int ,char ** argv)
//{

//    // dataset july
//    velodyne_rpy = vcg::Point3f(3.151592653589793,0.1323284641020683,0.);
//    left_rpy = vcg::Point3f(-1.5708, 0, -1.5708);

//    lidar2imu = rpy2mat(velodyne_rpy);
//    camera2imu = rpy2mat(left_rpy);

//    // dataset november
//    lidar2imu = vcg::Matrix44f(T_scan_imu);
//    camera2imu = vcg::Matrix44f(T_cam_imu);


//    scale = 1.0;


//    // read imu values
//    readImuTS(std::string(argv[1]));
//    // read image time stamps
//    readImagesTime(std::string(argv[2]));
//    // read images aligned cameras
//    import_from_bundle(argv[3]);

//    align_bins(std::string(argv[4]));
//    vcg::tri::io::ExporterPLY<MyMesh>::Save(::allLerpCams,"all_lerp_cams.ply",vcg::tri::io::Mask::IOM_ALL);
//    return 0;
//}

void readLidarTranslations(char * folder){
    int i=0;
    while(readBinTranslation((std::string(folder)+"/"+std::to_string(i)+".bin"))){
        i++;
    }
}

void bins2ply(char * folder){
    int i=0;
    while( readBinMeshNoSfm((std::string(folder)+"/"+std::to_string(i)+".bin"),tessellate)){
        i+=bin_step;
        converted_bins.push_back(std::to_string(i)+".bin");
    }
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

int create_mlp(char * folder, bool tessellate ){
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
    for(int i=0; i < converted_bins.size(); i++)
    of <<"<MLMesh label=\""<< converted_bins[i] <<"\" filename=\""<< "LIDAR/"<<converted_bins[i]<<mname <<"\">\n<MLMatrix44> 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 </MLMatrix44> </MLMesh>" << std::endl;
    of <<"</MeshGroup>\n";

    of << "<RasterGroup>" << std::endl;
    for(int i=0; i <  image_names.size();i++){
    p = positionFromBin(image_stamps[i]);
    p+=lidar2imu.GetColumn3(3);// to imu
    p-=camera2imu.GetColumn3(3);// to camera


    //vcg::Quaterniond qrot = rotationFromImu(image_stamps[i]);
    // qrot.ToMatrix(imuFrame);

    // imu2camera = vcg::Inverse(camera2imu);
    // imu2camera.SetColumn(3,vcg::Point4d(0,0,0,1));

    // cameraFrame = imuFrame*imu2camera;

     int il = -1;
     vcg::Quaterniond qrotL = rotationFromBin(image_stamps[i]+10.0,il);
     vcg::Matrix44d lidarFrame;
     qrotL.ToMatrix(lidarFrame);

     printf("lidar n. %d\n",il);


     imu2camera = vcg::Inverse(camera2imu);
     imu2camera.SetColumn(3,vcg::Point4d(0,0,0,1));
     lidar2imu.SetColumn(3,vcg::Point4d(0,0,0,1));




     cameraFrame = lidarFrame*lidar2imu*imu2camera;
     cameraFrame = cameraFrame*RR.SetRotateDeg(180,vcg::Point3d(1,0,0));

//     vcg::Point3d s = skew(camera2imu);
//     printf("camera2imu skew: %f %f %f\n",s[0],s[1],s[2]);

     // PATCH for Meshlab which goes nuts for up vector close to Z
     RR.SetRotateDeg(90,vcg::Point3d(1,0,0));
     cameraFrame = RR*cameraFrame;
     p = RR*p;


    // meshlab stores the camera as the rotation matrix and .
    // 3D to camera: rot*(point-p) (p camera center)
    // the camera frame is encoded in the rows of rot

   rot = cameraFrame;
   rot.transposeInPlace();
    of << "<MLRaster label=\""  << image_names[i].c_str() << "\">\n"
       << "<VCGCamera FocalMm=\""<< camera.focal << "\" RotationMatrix=\""<<
          rot[0][0] << " " << rot[0][1] <<" "<< rot[0][2] << " " << 0<< " "<<
          rot[1][0] << " " << rot[1][1] <<" "<< rot[1][2] << " " << 0<< " "<<
          rot[2][0] << " " << rot[2][1] <<" "<< rot[2][2] << " " << 0<< " "<<
          rot[3][0] << " " << rot[3][1] <<" "<< rot[3][2] << " " << 1<<
          "\" ViewportPx=\""<< camera.vp[0] << " " << camera.vp[1] <<"\" LensDistortion=\"0 0\"  BinaryData=\"0\"  CameraType=\"0\" TranslationVector=\""<<-p[0] << " " <<-p[1] << " "<<-p[2] <<" 1\" PixelSizeMm=\""<< 1.f   <<" " << 1.f <<"\""<<std::endl;
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

    // read image time stamps
    readImagesTime(std::string(argv[2]));

    // convert bin to plys
    bins2ply(argv[3]);

    // create MLP
    create_mlp(argv[3],tessellate);

    stats_lidar(argv[3]);
    return 0;
}
