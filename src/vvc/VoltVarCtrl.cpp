
#include "VoltVarCtrl.hpp"

#include <iostream>
#include <fstream>
//#include <chrono>
#include <ctime>
#include "CLogger.hpp"
#include "Messages.hpp"
#include "CTimings.hpp"
#include "CDeviceManager.hpp"
#include "CGlobalPeerList.hpp"
#include "gm/GroupManagement.hpp"
#include "CGlobalConfiguration.hpp"
//#include "device/COpenDssAdapter.hpp"

#include <sstream>
#include <string>
#include <vector>
#include <math.h>

#include <boost/foreach.hpp>
#include <boost/bind.hpp>
#include <boost/asio/error.hpp>
#include <boost/system/error_code.hpp>
#include <boost/range/adaptor/map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include <boost/foreach.hpp>
#include <boost/bind.hpp>
#include <boost/asio/error.hpp>
#include <boost/system/error_code.hpp>
#include <boost/range/adaptor/map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include <armadillo>

namespace freedm {
	
namespace broker {
		
namespace vvc { 

namespace {
/// This file's logger.
CLocalLogger Logger(__FILE__);
}

///////////////////////////////////////////////////////////////////////////////
/// VVCAgent
/// @description: Constructor for the VVC module
/// @pre: Posix Main should register read handler and invoke this module
/// @post: Object is initialized and ready to run 
/// @param uuid_: This object's uuid
/// @param broker: The Broker
/// @limitations: None
///////////////////////////////////////////////////////////////////////////////
VVCAgent::VVCAgent()
    : ROUND_TIME(boost::posix_time::milliseconds(CTimings::Get("LB_ROUND_TIME")))
    , REQUEST_TIMEOUT(boost::posix_time::milliseconds(CTimings::Get("LB_REQUEST_TIMEOUT")))
{

  Logger.Trace << __PRETTY_FUNCTION__ << std::endl;

  m_RoundTimer = CBroker::Instance().AllocateTimer("vvc");
  m_WaitTimer = CBroker::Instance().AllocateTimer("vvc");

	t1 = boost::posix_time::second_clock::local_time();
	std::cout << "DPF started at: " << t1 << std::endl;
	
	timer_overvolt = boost::posix_time::second_clock::local_time();
	timer_undervolt = boost::posix_time::second_clock::local_time();
	flag = 0;

	sysdata sysinfo = load_system_data();
  	Dl_old = sysinfo.Dl;
  	pload_vector.zeros(1014,1);
  	m_Synchronized = 0;
	
	PV_node << 7 << 20 << 23 << 28 << 30 << 36 << 44 << 55 << 60 << 73 << 77 << 102 << 110 << 114 << 116 << arma::endr;

}
VVCAgent::~VVCAgent()
{
}
			
////////////////////////////////////////////////////////////
/// Run
/// @description Main function which initiates the algorithm
/// @pre: Posix Main should invoke this function
/// @post: Triggers the vvc algorithm by calling VVCManage()
/// @limitations None
/////////////////////////////////////////////////////////
int VVCAgent::Run()
{
  std::cout<< " --------------------VVC ---------------------------------" << std::endl; 
  CBroker::Instance().Schedule("vvc",
      boost::bind(&VVCAgent::FirstRound, this, boost::system::error_code()));
  Logger.Info << "VVC is scheduled for the next phase." << std::endl;
  return 0;
}


///////////////////////////////////////////////////////////////////////////////
/// HandleIncomingMessage
/// "Downcasts" incoming messages into a specific message type, and passes the
/// message to an appropriate handler.
/// @pre None
/// @post The message is handled by the target handler or a warning is
///     produced.
/// @param m the incoming message
/// @param peer the node that sent this message (could be this DGI)
///////////////////////////////////////////////////////////////////////////////
void VVCAgent::HandleIncomingMessage(boost::shared_ptr<const ModuleMessage> m, CPeerNode peer)
{
    if(m->has_volt_var_message())
    {
        VoltVarMessage vvm = m->volt_var_message();
        if(vvm.has_voltage_delta_message())
        {
            HandleVoltageDelta(vvm.voltage_delta_message(), peer);
        }
        else if(vvm.has_line_readings_message())
        {
            HandleLineReadings(vvm.line_readings_message(), peer);
        }
	else if(vvm.has_gradient_message())
        {
            HandleGradient(vvm.gradient_message(), peer);
        }
        else
        {
            Logger.Warn << "Dropped unexpected volt var message: \n" << m->DebugString();
        }
    }
    else if(m->has_group_management_message())
    {
        gm::GroupManagementMessage gmm = m->group_management_message();
        if(gmm.has_peer_list_message())
        {
            HandlePeerList(gmm.peer_list_message(), peer);
        }
        else
        {
            Logger.Warn << "Dropped unexpected group management message:\n" << m->DebugString();
        }
    }
    else
    {
        Logger.Warn<<"Dropped message of unexpected type:\n" << m->DebugString();
    }
}

void VVCAgent::HandleVoltageDelta(const VoltageDeltaMessage & m, CPeerNode peer)
{
    Logger.Trace << __PRETTY_FUNCTION__ << std::endl;
    Logger.Notice << "Got VoltageDelta from: " << peer.GetUUID() << std::endl;
    Logger.Notice << "CF "<<m.control_factor()<<" Phase "<<m.phase_measurement()<<std::endl;
}

void VVCAgent::HandleLineReadings(const LineReadingsMessage & m, CPeerNode peer)
{
    Logger.Trace << __PRETTY_FUNCTION__ << std::endl;
    Logger.Notice << "Got Line Readings from "<< peer.GetUUID() << std::endl;
}

void VVCAgent::HandleGradient(const GradientMessage & m, CPeerNode peer)
{
    Logger.Trace << __PRETTY_FUNCTION__ << std::endl;
    Logger.Notice << "Got Gradients from "<< peer.GetUUID() << std::endl;
    Logger.Notice << "size of vector "<< m.gradient_value_size() << std::endl;
    Logger.Notice << "the 1st element = " << m.gradient_value(0) <<std::endl;

	/**** NEW Code ****/
	if (peer.GetUUID() == "arrogate:5002")
	{
		for (int i = 0; i < m.gradient_value_size(); i++)
		pload_vector(i,0) = m.gradient_value(i);
		std::cout << "Loaded Pload from: " << peer.GetUUID() << std::endl;
		m_Synchronized = 1;
	}
    
}


// HandlePeerlist Implementation
void VVCAgent::HandlePeerList(const gm::PeerListMessage & m, CPeerNode peer)
{
	Logger.Trace << __PRETTY_FUNCTION__ << std::endl;
	Logger.Notice << "Updated Peer List Received from: " << peer.GetUUID() << std::endl;
	
	m_peers.clear();
	
	m_peers = gm::GMAgent::ProcessPeerList(m);
	m_leader = peer.GetUUID();
}

// Preparing Messages
ModuleMessage VVCAgent::VoltageDelta(unsigned int cf, float pm, std::string loc)
{
	VoltVarMessage vvm;
	VoltageDeltaMessage *vdm = vvm.mutable_voltage_delta_message();
	vdm -> set_control_factor(cf);
	vdm -> set_phase_measurement(pm);
	vdm -> set_reading_location(loc);
	return PrepareForSending(vvm,"vvc");
}

ModuleMessage VVCAgent::LineReadings(std::vector<float> vals)
{
	VoltVarMessage vvm;
	std::vector<float>::iterator it;
	LineReadingsMessage *lrm = vvm.mutable_line_readings_message();
	for (it = vals.begin(); it != vals.end(); it++)
	{
		lrm -> add_measurement(*it);
	}
	lrm->set_capture_time(boost::posix_time::to_simple_string(boost::posix_time::microsec_clock::universal_time()));
	return PrepareForSending(vvm,"vvc");
}

ModuleMessage VVCAgent::Gradient(arma::mat grad)
{
	VoltVarMessage vvm;
	unsigned int idx;
        GradientMessage *grdm = vvm.mutable_gradient_message();
	for (idx = 0; idx < grad.n_rows; idx++)
	{
		grdm -> add_gradient_value(grad(idx));
	}
	grdm->set_gradient_capture_time(boost::posix_time::to_simple_string(boost::posix_time::microsec_clock::universal_time()));
	return PrepareForSending(vvm,"vvc");
}


ModuleMessage VVCAgent::PrepareForSending(const VoltVarMessage& message, std::string recipient)
{
    Logger.Trace << __PRETTY_FUNCTION__ << std::endl;
    ModuleMessage mm;
    mm.mutable_volt_var_message()->CopyFrom(message);
    mm.set_recipient_module(recipient);
    return mm;
}
// End of Preparing Messages

///////////////////////////////////////////////////////////////////////////////
/// FirstRound
/// @description The code that is executed as part of the first VVC
///     each round.
/// @pre None
/// @post if the timer wasn't cancelled this function calls the first VVC.
/// @param error The reason this function was called.
///////////////////////////////////////////////////////////////////////////////
void VVCAgent::FirstRound(const boost::system::error_code& err)
{
  Logger.Trace << __PRETTY_FUNCTION__ << std::endl;
  
  if(!err)
  {
    CBroker::Instance().Schedule("vvc", 
	boost::bind(&VVCAgent::VVCManage, this, boost::system::error_code()));
  }
  else if(err == boost::asio::error::operation_aborted)
  {
    Logger.Notice << "VVCManage Aborted" << std::endl;
  }
  else
  {
    Logger.Error << err << std::endl;
    throw boost::system::system_error(err);
  }
			
}


////////////////////////////////////////////////////////////
/// VVCManage
/// @description: Manages the execution of the VVC algorithm 
/// @pre: 
/// @post: 
/// @peers 
/// @limitations
/////////////////////////////////////////////////////////
void VVCAgent::VVCManage(const boost::system::error_code& err)
{
    Logger.Trace << __PRETTY_FUNCTION__ << std::endl;

    if(!err)
    {
      ScheduleNextRound();
      ReadDevices();
      vvc_main();
     
      
    }
        
    else if(err == boost::asio::error::operation_aborted)
    {
        Logger.Notice << "VVCManage Aborted" << std::endl;
    }
    else
    {
        Logger.Error << err << std::endl;
        throw boost::system::system_error(err);
    }
}

///////////////////////////////////////////////////////////////////////////////
/// ScheduleNextRound
/// @description Computes how much time is remaining and if there isn't enough
///     requests the VVC that will run next round.
/// @pre None
/// @post VVCManage is scheduled for this round OR FirstRound is scheduled
///     for next time.
///////////////////////////////////////////////////////////////////////////////
void VVCAgent::ScheduleNextRound()
{
    Logger.Trace << __PRETTY_FUNCTION__ << std::endl;

    if(CBroker::Instance().TimeRemaining() > ROUND_TIME + ROUND_TIME)
    {
        CBroker::Instance().Schedule(m_RoundTimer, ROUND_TIME,
            boost::bind(&VVCAgent::VVCManage, this, boost::asio::placeholders::error));
        Logger.Info << "VVCManage scheduled in " << ROUND_TIME << " ms." << std::endl;
    }
    else
    {
        CBroker::Instance().Schedule(m_RoundTimer, boost::posix_time::not_a_date_time,
            boost::bind(&VVCAgent::FirstRound, this, boost::asio::placeholders::error));
        Logger.Info << "VVCManage scheduled for the next phase." << std::endl;
    }
}


///////////////////////////////////////////////////////////////////////////////
/// ReadDevices
/// @description Reads the device state and updates the appropriate member vars.
/// @pre None
/// @post 
///////////////////////////////////////////////////////////////////////////////
void VVCAgent::ReadDevices()
{
    Logger.Trace << __PRETTY_FUNCTION__ << std::endl;

    float generation = device::CDeviceManager::Instance().GetNetValue("Drer", "generation");
    float storage = device::CDeviceManager::Instance().GetNetValue("Desd", "storage");
    float load = device::CDeviceManager::Instance().GetNetValue("Load", "drain");

    m_Gateway = device::CDeviceManager::Instance().GetNetValue("Sst", "gateway");
    m_NetGeneration = generation + storage - load;
}



///////////////////////////////////////////////////////////////////////////////
/// Main VVO Code with 4 VR module
/// @description Reads the device states and provides Qinj & VR tap commands.
/// @pre None
/// @post 
///////////////////////////////////////////////////////////////////////////////



void VVCAgent::vvc_main()
{
using namespace arma;
using namespace std;
	
//Prepare para for DPF
sysdata sysinfo = load_system_data();
int Ldl = sysinfo.Dl.n_rows;
int Wdl = sysinfo.Dl.n_cols;
double v0 = sysinfo.vo;
cout << "Dl dimension:"<< Ldl <<"*"<< Wdl << endl;//Matrix Dl in Matlab

//mat Dl = sysinfo.Dl;
cx_mat Z = sysinfo.Z;
cout << "Z dimensions:" << Z.n_rows << "*" << Z.n_cols << endl;
mat Dl = Dl_old;
mat Dl_pv = sysinfo.Dl_pv;

//mat lcd_vr1, lcd_vr2, lcd_vr3, lcd_vr4;

//lcd_vr1 << 1 << 1 << 1 << endr;
//lcd_vr2 << 1 << 1 << 1 << endr;
//lcd_vr3 << 1 << 0 << 0 << endr;
//lcd_vr4 << 1 << 0 << 1 << endr;

//document the original node number in sequence
int j, ja, jb, jc, ia, ib, ic;

j = 1;
ja = 0;
jb = 0;
jc = 0;

int cnt_nodes = 0;
int Lla = 0;
int Llb = 0;
int Llc = 0;
for (int i = 0; i < Ldl; ++i)
{
  if ((int)Dl(i, 0) != 0)
  cnt_nodes = cnt_nodes + 1;
  if ((int)Dl(i, 6) != 0)
    Lla = Lla + 1;
  if ((int)Dl(i, 8) != 0)
    Llb = Llb + 1;
  if ((int)Dl(i, 10) != 0)
    Llc = Llc + 1;
}
cnt_nodes = cnt_nodes + 1;//No of nodes= No of branches +1

/**** NEW Code ****/
arma::mat Pload_mat;
if (m_Synchronized)
{
	mat vt1, vt2, vt3, vt4;
	Pload_mat = pload_vector;
	cout << "size of Pload_vec: " << pload_vector.n_rows << endl;
	cout << "size of Pload_mat: " << Pload_mat.n_rows << endl;
	Dl.col(6) = Pload_mat.rows(0,166);
	Dl.col(8) = Pload_mat.rows(167,333);
	Dl.col(10) = Pload_mat.rows(334,500);
	Dl.col(7) = Pload_mat.rows(501, 667);
	Dl.col(9) = Pload_mat.rows(668, 834);
	Dl.col(11) = Pload_mat.rows(835, 1001);
	vt1 = Pload_mat.rows(1002, 1004);
	vt2 = Pload_mat.rows(1005, 1007);
	vt3 = Pload_mat.rows(1008, 1010);
	vt4 = Pload_mat.rows(1011, 1013);
	Dl(0, span(13, 15)) = vt1.t();
	Dl(12, span(13, 15)) = vt2.t();
	Dl(34, span(13, 15)) = vt3.t();
	Dl(153, span(13, 15)) = vt4.t(); 
	cout << "Dl size: " << Dl.n_rows << endl;
	m_Synchronized = 0;
	flag1 = 1;
}
else
{
	cout << "Couldn't load Pload_vector" << endl;
	Dl.cols(6,11) = Dl_old.cols(6,11);
	Dl(0, span(13, 15)) = Dl_old(0, span(13, 15));
	Dl(12, span(13, 15)) = Dl_old(12, span(13, 15));
	Dl(34, span(13, 15)) = Dl_old(34, span(13, 15));
	Dl(153, span(13, 15)) = Dl_old(153, span(13, 15));
}
	
/**** till here ****/

int Nnum = cnt_nodes;
mat aVR1, aVR2, aVR3, aVR4;
aVR1 = Dl(0, span(13, 15));
aVR2 = Dl(12, span(13, 15));
aVR3 = Dl(34, span(13, 15));
aVR4 = Dl(153, span(13, 15));

//int timer_overvolt = 0;
//int timer_undervolt = 0;
mat Dl_vr;
mat Qcmd;
double Qmax_his, Qmin_his;

//auto start = std::chrono::system_clock::now();
//VCtrig vctrig = VoltReg_trigger(Dl, Z, Nnum, timer_overvolt, timer_undervolt);
VCtrig vctrig = VoltReg_trigger(Dl, Z, Nnum, timer_overvolt, timer_undervolt, flag);
flag = vctrig.flag;
if (flag == 1)
{
	timer_overvolt = vctrig.timer_overv;
	timer_undervolt = vctrig.timer_underv;
}
else
{
	timer_overvolt = boost::posix_time::second_clock::local_time();
	timer_undervolt = boost::posix_time::second_clock::local_time();
}
cout << "Voltage module trig = " << vctrig.trigger << endl;
if (vctrig.trigger == 1)
{
	VRreg vrreg = volt_reg_solver(Dl, Nnum, v0);
	Dl_vr = vrreg.Dl_vr;
	mat re_vr = vrreg.re_vr;
	aVR1 = Dl_vr(0, span(13, 15));
	aVR2 = Dl_vr(12, span(13, 15));
	aVR3 = Dl_vr(34, span(13, 15));
	aVR4 = Dl_vr(153, span(13, 15));
	//cout << "aVR1: " << aVR1 << endl;
	//cout << "aVR2: " << aVR2 << endl;
	//cout << "aVR3: " << aVR3 << endl;
	//cout << "aVR4: " << aVR4 << endl;
	cout << "Voltage Regulation (4 VRs) Completed! \n";
}
else
{
	cout << "Voltage Regulation (4 VRs) Skipped \n";
	Dl_vr = Dl;
	aVR1 = Dl_vr(0, span(13, 15));
	aVR2 = Dl_vr(12, span(13, 15));
	aVR3 = Dl_vr(34, span(13, 15));
	aVR4 = Dl_vr(153, span(13, 15));
}
double bkva = sysinfo.bkva;
gradV gradv = grad_Volt_opt_solve(Dl_vr, v0, Z, bkva, Dl_pv, PV_node);
//cout << "Dl_vr: \n" << Dl_vr(span::all, span(13,15)) << endl;
mat Dl_volt = gradv.Dl_volt;

int trigger1 = gradv.trigger1;

if (trigger1 == 1)
{
	gradVAR gradq = var_opt_loss_solver(Dl_volt, v0, Z, bkva, Dl_pv, PV_node);
	Qcmd = gradq.Qcmd;
	cout << "trigger1 == 1" << endl;
	//cout << "Qcmd: \n" << Qcmd << endl;
}
else
{
	gradVAR gradq = var_opt_loss_solver(Dl_vr, v0, Z, bkva, Dl_pv, PV_node);
	Qcmd = gradq.Qcmd;
	cout << "trigger1 == 0" << endl;
}

//cout << "Qcmd: \n" << Qcmd << endl;
//cout << "size of Qcmd: " << Qcmd.n_rows << "*" << Qcmd.n_cols << endl;
uvec temp = find(Qcmd != 0);
if (temp.is_empty())
{
	Qmax_his = 0;
	Qmin_his = 0;
}
else
{
	//Qmax_his = max(max(Qcmd(find()))
	mat tempA = Qcmd.col(0);
	mat tempB = Qcmd.col(1);
	mat tempC = Qcmd.col(2);
	uvec idx_A = find(tempA != 0);
	uvec idx_B = find(tempB != 0);
	uvec idx_C = find(tempC != 0);
	mat Qmax_his_abc, Qmin_his_abc;
	Qmax_his_abc << max(tempA(idx_A)) << max(tempB(idx_B)) << max(tempC(idx_C)) << endr;
	Qmin_his_abc << min(tempA(idx_A)) << min(tempB(idx_B)) << min(tempC(idx_C)) << endr;
	Qmax_his = max(max(Qmax_his_abc));
	Qmin_his = min(min(Qmin_his_abc));
}
//cout << "Qmax_his: " << Qmax_his << "\nQmin_his: " << Qmin_his << endl;

//mat Dl_check = Dl_vr;
Dl_vr.col(7) = Qcmd.col(0);
Dl_vr.col(9) = Qcmd.col(1);
Dl_vr.col(11) = Qcmd.col(2);

//cout << "Qcmd: \n" << Qcmd << endl;
//cout << "Reached almost the end!" << endl;

cout << "Qcmd rows: " << Qcmd.n_rows << endl;
VPQ dpf_re = DPF_return123(Dl_vr, Z);
mat Vpolar = dpf_re.Vpolar;

mat Vresult, Vpolar_a, Vpolar_b, Vpolar_c, Vpolar_tempa, Vpolar_tempb, Vpolar_tempc;
Vpolar_tempa = Vpolar.col(0);
Vpolar_a = Vpolar_tempa.elem(find(Vpolar.col(0) != 0));
Vpolar_tempb = Vpolar.col(2);
Vpolar_b = Vpolar_tempb.elem(find(Vpolar.col(2) != 0));
Vpolar_tempc = Vpolar.col(4);
Vpolar_c = Vpolar_tempc.elem(find(Vpolar.col(4) != 0));
//cout << "Vpolar: \n" << Vpolar << endl;
mat VV = join_cols(Vpolar_a, join_cols(Vpolar_b, Vpolar_c));
mat V_var = var(VV);
//cout << "\nV_variance: " << (V_var*10000) << endl;
//cout << "V: " << VV << endl;
//system("pause");


if (flag1)
{
		BOOST_FOREACH(CPeerNode peer, m_peers | boost::adaptors::map_values)
        	{		
			
		    	mat S2, temp1, temp2, temp3, temp4;
			temp1 = aVR1.t();
			temp2 = aVR2.t();
			temp3 = aVR3.t();
			temp4 = aVR4.t();
			S2 << Dl_vr(19,7) << Dl_vr(22,7) << Dl_vr(90,11) << Dl_vr(85,9) << Dl_vr(138,7) << Dl_vr(145,11) << Dl_vr(49,7) << Dl_vr(49,9) << Dl_vr(49,11) << Dl_vr(128,7) << Dl_vr(128,9) << Dl_vr(128,11) << Dl_vr(63,7) << Dl_vr(63,9) << Dl_vr(63,11) << Dl_vr(70,7) << Dl_vr(70,9) << Dl_vr(70,11) << Dl_vr(6,7) << Dl_vr(40,11) << Dl_vr(149,9) << Dl_vr(33,7) << Dl_vr(29,11) << endr;
			
			S2 = S2.t();
			//S2 = join_cols(join_cols(Qcmd.col(0), Qcmd.col(1)), Qcmd.col(2));
			S2 = join_cols(S2, temp1);
			S2 = join_cols(S2, temp2);
			S2 = join_cols(S2, temp3);
			S2 = join_cols(S2, temp4);
			cout << "S2: \n" << S2 << endl;
			cout << "Size of Gradient (S2) sent to slaves is: " << S2.n_rows << endl;
			if (peer.GetUUID() != "arrogate:5002" && peer.GetUUID() != GetUUID())
			{

			    	//cout << "Size of S2: " << S2.n_rows << endl;
			    	cout << "Sending to: " << peer.GetUUID() << " Gradient Data!" << endl;
				ModuleMessage mg = Gradient(S2);
			    	peer.Send(mg);
				flag1 = 0;
			}
	    
		}
}

/***************** Testing Sending Capability of Opendss Adapter *********************/
/*mat S2;
//S2 << Dl_osize(1,7) << Dl_osize(2,7) << Dl_osize(3,7) << Dl_osize(4,7) << Dl_osize(6,7) << Dl_osize(7,7) << Dl_osize(8,7) << Dl_osize(1,9) << Dl_osize(2,9) << Dl_osize(3,9) << Dl_osize(4,9) << Dl_osize(6,9) << Dl_osize(7,9) << Dl_osize(8,9) << Dl_osize(1,11)<< Dl_osize(2,11)<< Dl_osize(3,11)<< Dl_osize(4,11)<< Dl_osize(6,11)<< Dl_osize(7,11)<< Dl_osize(8,11)<< endr;

for (int i = 0; i < Qcmd.n_rows; i++)
	for (int j = 0; j < Qcmd.n_cols; j++)
		S2 << Qcmd(i,j);

S2 << aVR1(0,0) << aVR1(0,1) << aVR1(0,2) << aVR2(0,0) << aVR2(0,1) << aVR2(0,2) << aVR3(0,0) << aVR3(0,1) << aVR3(0,2) << aVR4(0,0) << aVR4(0,1) << aVR4(0,2) << endr;
	
//S2 << Dl_osize(1,7) << Dl_osize(1,9) << Dl_osize(1,11) << Dl_osize(2,7) << Dl_osize(2,9) << Dl_osize(2,11) << Dl_osize(3,7) << Dl_osize(3,9) << Dl_osize(3,11) << Dl_osize(4,7) << Dl_osize(4,9) << Dl_osize(4,11) << Dl_osize(6,7) << Dl_osize(6,9) << Dl_osize(6,11) << Dl_osize(7,7) << Dl_osize(7,9) << Dl_osize(7,11) << Dl_osize(8,7) << Dl_osize(8,9) << Dl_osize(8,11) << endr;
S2 = S2.t();

//Code added by Valli on 6/19/18
//Conversion to string to be sent to opendss
ostringstream stream;

for (int i = 0; i < S2.n_rows; i++)
{
	stream << S2(i,0) << '\n';
}

//cout << "S2: \n" << S2 << endl;
cout << "Stream: \n" << stream.str() << endl;
std::string command = stream.str(); // generic command should be changed
device::COpenDssAdapter::sendCommand(command);    //test sendop*/

/***************** Testing Sending Capability of Opendss Adapter ENDS HERE*********************/

//Dl_check
//save three phase branch data separately
//cx_mat brn_a = cx_mat(zeros(Y_return.Lnum_a, Y_return.brnches.n_cols), zeros(Y_return.Lnum_a, Y_return.brnches.n_cols));
//cx_mat brn_b = cx_mat(zeros(Y_return.Lnum_b, Y_return.brnches.n_cols), zeros(Y_return.Lnum_b, Y_return.brnches.n_cols));
//cx_mat brn_c = cx_mat(zeros(Y_return.Lnum_c, Y_return.brnches.n_cols), zeros(Y_return.Lnum_c, Y_return.brnches.n_cols));
//				
//ja = 0;
//jb = 0;
//jc = 0;
//
//for (int i = 0; i < Lbr && ja<Y_return.Lnum_a && jb<Y_return.Lnum_b && jc<Y_return.Lnum_c; ++i)
//{
//  if (abs(Y_return.brnches(i, 2)) != 0)
//  {
//    brn_a.row(ja) = Y_return.brnches.row(i);
//    ++ja;
//  }
//  if (abs(Y_return.brnches(i, 3)) != 0)
//  {
//    brn_b.row(jb) = Y_return.brnches.row(i);
//    ++jb;
//  }
//  if (abs(Y_return.brnches(i, 4)) != 0)
//  {
//    brn_c.row(jc) = Y_return.brnches.row(i);
//    ++jc;
//  }
//}

 
//  // send messages to slaves
//if (Ploss_osize < Ploss_orig)// grad message will NOT be sent to slaves if loss is not reduced
//{
//  BOOST_FOREACH(CPeerNode peer, m_peers | boost::adaptors::map_values)
//        {
//      	    
//	    ModuleMessage mm = VoltageDelta(2, 3.0, "NCSU");
//            peer.Send(mm);
//	    
//	   
//	    mat S2;
//	    S2 << Dl(1,7) << Dl(2,7) << Dl(3,7) << Dl(4,7) << Dl(6,7) << Dl(7,7) << Dl(8,7) << Dl(1,9) << Dl(2,9) << Dl(3,9) << Dl(4,9) << Dl(6,9) << Dl(7,9) << Dl(8,9) << Dl(1,11)<< Dl(2,11)<< Dl(3,11)<< Dl(4,11)<< Dl(6,11)<< Dl(7,11)<< Dl(8,11)<< endr;
//	    S2 = S2.t();
//	    ModuleMessage mg = Gradient(S2);
//	    peer.Send(mg);
//	    
//	}
//}
//else
//{
//}
//  
//  break;

//   send messages to slaves
//if (Ploss_osize < Ploss_orig)// grad message will NOT be sent to slaves if loss is not reduced
//{
//  BOOST_FOREACH(CPeerNode peer, m_peers | boost::adaptors::map_values)
//        {
//      	    
//	    ModuleMessage mm = VoltageDelta(2, 3.0, "Gradients reversed!");
//            peer.Send(mm);
//	    
//	   
//	    mat S2;
//	    S2 << Dl(1,7) << Dl(2,7) << Dl(3,7) << Dl(4,7) << Dl(6,7) << Dl(7,7) << Dl(8,7) << Dl(1,9) << Dl(2,9) << Dl(3,9) << Dl(4,9) << Dl(6,9) << Dl(7,9) << Dl(8,9) << Dl(1,11)<< Dl(2,11)<< Dl(3,11)<< Dl(4,11)<< Dl(6,11)<< Dl(7,11)<< Dl(8,11)<< endr;
//	    S2 = S2.t();
//	    ModuleMessage mg = Gradient(S2);
//	    peer.Send(mg);
//	}
//}
//else
//{
//}
//	
//  break;

//Time duration calculation
t2 = boost::posix_time::second_clock::local_time();

//t2(boost::posix_time::second_clock::local_time());
td = t2 - t1;
cout << "Time elapsed is: " << td.seconds() << endl;
t1 = t2;

Dl_old = Dl;

}// end of vvc_main()


}//namespace vvc
}// namespace broker
}// namespace freedm
//}
