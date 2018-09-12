
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

	i = 0; //The load iterator is initiated to 0 once at beginning
	std::ifstream file1, file2;
	file1.open("r_load.dat");
	load_profile.load(file1);
	file1.close();
	file2.open("r_pv.dat");
	pv_profile.load(file2);
	file2.close();
	std::cout << "size of load_profile: " << load_profile.n_rows << "*" << load_profile.n_cols << std::endl;
	std::cout << "size of pv_profile: " << pv_profile.n_rows << "*" << pv_profile.n_cols << std::endl;
	t0 = boost::posix_time::microsec_clock::local_time();
	t1 = boost::posix_time::second_clock::local_time();
	std::cout << "DPF started at: " << t1 << std::endl;

	l = 0;
	l1 = 0;
	
	//grad_slav.zeros(513,1);

	grad_slav1.zeros(3,1);
	grad_slav2.zeros(3,1);
	grad_slav3.zeros(3,1);
	grad_slav4.zeros(3,1);
	grad_slav5.zeros(3,1);
	grad_slav6.zeros(3,1);
	grad_slav7.zeros(3,1);
	grad_slav8.zeros(3,1);
	grad_slav9.zeros(3,1);
	grad_slav10.zeros(3,1);
	grad_slav11.zeros(3,1);
	grad_slav12.zeros(3,1);
	grad_slavvr1.ones(3,1);
	grad_slavvr2.ones(3,1);
	grad_slavvr3.ones(3,1);
	grad_slavvr4.ones(3,1);
	
	m_Synchronized1 = 0;
	m_Synchronized2 = 0;
	m_Synchronized3 = 0;
	m_Synchronized4 = 0;
	m_Synchronized5 = 0;
	m_Synchronized6 = 0;
	m_Synchronized7 = 0;
	m_Synchronized8 = 0;
	m_Synchronized9 = 0;
	m_Synchronized10 = 0;
	m_Synchronized11 = 0;
	m_Synchronized12 = 0;
	m_Synchronizedvr1 = 0;
	m_Synchronizedvr2 = 0;
	m_Synchronizedvr3 = 0;
	m_Synchronizedvr4 = 0;
	
	PV_node << 7 << 20 << 23 << 28 << 30 << 36 << 44 << 55 << 60 << 73 << 77 << 102 << 110 << 114 << 116 << arma::endr;
	
	m_Synchronized = 0;

	time = arma::zeros<arma::vec>(60000);
	min_v = arma::zeros<arma::vec>(300);
	max_v = arma::zeros<arma::vec>(300);
	total_loss = arma::zeros<arma::vec>(300);
	delay = arma::zeros<arma::vec>(60000);

	m_S1 = arma::zeros<arma::vec>(600);
	m_S2 = arma::zeros<arma::vec>(600);
	m_S3 = arma::zeros<arma::vec>(600);
	m_S4 = arma::zeros<arma::vec>(600);
	m_S5 = arma::zeros<arma::vec>(600);
	m_S6 = arma::zeros<arma::vec>(600);
	lx = 0;

	sysdata sysinfo = load_system_data();
	Dl_old = sysinfo.Dl;
	Dl_pv_old = sysinfo.Dl_pv;
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
	if (peer.GetUUID() == "arrogate:5003")
	{
		for (int i = 0; i < m.gradient_value_size(); i++)
		grad_slav1(i,0) = m.gradient_value(i);
		std::cout << "Loaded Pload from: " << peer.GetUUID() << std::endl;
		m_Synchronized1 = 1;
	}
	
	if (peer.GetUUID() == "arrogate:5004")
	{
		for (int i = 0; i < m.gradient_value_size(); i++)
		grad_slav2(i,0) = m.gradient_value(i);
		std::cout << "Loaded Pload from: " << peer.GetUUID() << std::endl;
		m_Synchronized2 = 1;
	}
	
	if (peer.GetUUID() == "arrogate:5005")
	{
		for (int i = 0; i < m.gradient_value_size(); i++)
		grad_slav3(i,0) = m.gradient_value(i);
		std::cout << "Loaded Pload from: " << peer.GetUUID() << std::endl;
		m_Synchronized3 = 1;
	}
	
	if (peer.GetUUID() == "arrogate:5006")
	{
		for (int i = 0; i < m.gradient_value_size(); i++)
		grad_slav4(i,0) = m.gradient_value(i);
		std::cout << "Loaded Pload from: " << peer.GetUUID() << std::endl;
		m_Synchronized4 = 1;
	}

	if (peer.GetUUID() == "arrogate:5007")
	{
		for (int i = 0; i < m.gradient_value_size(); i++)
		grad_slav5(i,0) = m.gradient_value(i);
		std::cout << "Loaded Pload from: " << peer.GetUUID() << std::endl;
		m_Synchronized5 = 1;
	}

	if (peer.GetUUID() == "arrogate:5008")
	{
		for (int i = 0; i < m.gradient_value_size(); i++)
		grad_slav6(i,0) = m.gradient_value(i);
		std::cout << "Loaded Pload from: " << peer.GetUUID() << std::endl;
		m_Synchronized6 = 1;
	}

	if (peer.GetUUID() == "arrogate:5009")
	{
		for (int i = 0; i < m.gradient_value_size(); i++)

		grad_slav7(i,0) = m.gradient_value(i);
		std::cout << "Loaded Pload from: " << peer.GetUUID() << std::endl;
		m_Synchronized7 = 1;
	}

	if (peer.GetUUID() == "arrogate:5010")
	{
		for (int i = 0; i < m.gradient_value_size(); i++)

		grad_slav8(i,0) = m.gradient_value(i);
		std::cout << "Loaded Pload from: " << peer.GetUUID() << std::endl;
		m_Synchronized8 = 1;
	}

	if (peer.GetUUID() == "arrogate:5011")
	{
		for (int i = 0; i < m.gradient_value_size(); i++)

		grad_slav9(i,0) = m.gradient_value(i);
		std::cout << "Loaded Pload from: " << peer.GetUUID() << std::endl;
		m_Synchronized9 = 1;
	}

	if (peer.GetUUID() == "arrogate:5012")
	{
		for (int i = 0; i < m.gradient_value_size(); i++)

		grad_slav10(i,0) = m.gradient_value(i);
		std::cout << "Loaded Pload from: " << peer.GetUUID() << std::endl;
		m_Synchronized10 = 1;
	}

	if (peer.GetUUID() == "arrogate:5013")
	{
		for (int i = 0; i < m.gradient_value_size(); i++)

		grad_slav11(i,0) = m.gradient_value(i);
		std::cout << "Loaded Pload from: " << peer.GetUUID() << std::endl;
		m_Synchronized11 = 1;
	}

	if (peer.GetUUID() == "arrogate:5014")
	{
		for (int i = 0; i < m.gradient_value_size(); i++)

		grad_slav12(i,0) = m.gradient_value(i);
		std::cout << "Loaded Pload from: " << peer.GetUUID() << std::endl;
		m_Synchronized12 = 1;
	}

	if (peer.GetUUID() == "arrogate:5015")
	{
		for (int i = 0; i < m.gradient_value_size(); i++)

		grad_slavvr1(i,0) = m.gradient_value(i);
		std::cout << "Loaded Pload from: " << peer.GetUUID() << std::endl;
		m_Synchronizedvr1 = 1;
	}

	if (peer.GetUUID() == "arrogate:5016")
	{
		for (int i = 0; i < m.gradient_value_size(); i++)
		grad_slavvr2(i,0) = m.gradient_value(i);
		std::cout << "Loaded Pload from: " << peer.GetUUID() << std::endl;
		m_Synchronizedvr2 = 1;
	}

	if (peer.GetUUID() == "arrogate:5017")
	{
		for (int i = 0; i < m.gradient_value_size(); i++)
		grad_slavvr3(i,0) = m.gradient_value(i);
		std::cout << "Loaded Pload from: " << peer.GetUUID() << std::endl;
		m_Synchronizedvr3 = 1;
	}

	if (peer.GetUUID() == "arrogate:5018")
	{
		for (int i = 0; i < m.gradient_value_size(); i++)
		grad_slavvr4(i,0) = m.gradient_value(i);
		std::cout << "Loaded Pload from: " << peer.GetUUID() << std::endl;
		m_Synchronizedvr4 = 1;
	}

	/**** till here ****/
    
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

cout << "Value of i is: " << i << endl;
cout << "Size of load_profile: " << load_profile.n_rows << endl;
cout << "Size of pv_profile: " << pv_profile.n_rows << endl;
cout << "Current load Multiplier: " << load_profile(i) << endl;
cout << "Current pv Multiplier: " << pv_profile(i) << endl;	
//Prepare para for DPF
sysdata sysinfo = load_system_data();
int Ldl = sysinfo.Dl.n_rows;
int Wdl = sysinfo.Dl.n_cols;
double v0 = sysinfo.vo;
cout << "Dl dimension:"<< Ldl <<"*"<< Wdl << endl;//Matrix Dl in Matlab

mat Dl_pv = sysinfo.Dl_pv;
mat Dl = sysinfo.Dl;
mat Dl_ld = sysinfo.Dl;
cx_mat Z = sysinfo.Z;
cout << "Z dimensions:" << Z.n_rows << "*" << Z.n_cols << endl;

y_re Y_return = form_Y_abc(Dl, sysinfo.Z, sysinfo.bkva, sysinfo.bkv);

//document the original node number in sequence
int j = 1;
int cnt_nodes = 0;

for (int k = 0; k < Ldl; ++k)
{
  if ((int)Dl(k, 0) != 0)
  cnt_nodes = cnt_nodes + 1;
}
cnt_nodes = cnt_nodes + 1;//No of nodes= No of branches +1
mat Node_f = zeros(1,cnt_nodes);

for (int k = 0; k < Ldl && j<cnt_nodes; ++k)
{
  if (	(int)Dl(k, 2) != 0  )
  {
    Node_f(0, j) = Dl(k, 2);
    ++j;
  }
}
Node_f = Node_f.st();

if (i==0)
{
	Dl.col(6) = load_profile(i)*Dl_ld.col(6) - pv_profile(i)*Dl_pv.col(6);
	Dl.col(8) = load_profile(i)*Dl_ld.col(8) - pv_profile(i)*Dl_pv.col(8);
	Dl.col(10) = load_profile(i)*Dl_ld.col(10) - pv_profile(i)*Dl_pv.col(10);
	Dl.col(7) = load_profile(i)*(Dl_ld.col(7) - Dl_pv.col(7));
	Dl.col(9) = load_profile(i)*(Dl_ld.col(7) - Dl_pv.col(7));
	Dl.col(11) = load_profile(i)*(Dl_ld.col(7) - Dl_pv.col(7));
}
else
{
	Dl.col(6) = load_profile(i)*Dl_ld.col(6) - pv_profile(i)*Dl_pv.col(6);
	Dl.col(8) = load_profile(i)*Dl_ld.col(8) - pv_profile(i)*Dl_pv.col(8);
	Dl.col(10) = load_profile(i)*Dl_ld.col(10) - pv_profile(i)*Dl_pv.col(10);
	Dl.col(7) = load_profile(i)*(Dl_ld.col(7) - Dl_pv.col(7));
	Dl.col(9) = load_profile(i)*(Dl_ld.col(7) - Dl_pv.col(7));
	Dl.col(11) = load_profile(i)*(Dl_ld.col(7) - Dl_pv.col(7));
	for (int z = 0; z < 15; z++)
	{
		Dl(PV_node(i), 7) = Dl_old(PV_node(i), 7);
		Dl(PV_node(i), 9) = Dl_old(PV_node(i), 9);
		Dl(PV_node(i), 11) = Dl_old(PV_node(i), 11);
	}
}

/*if (m_Synchronized)
{
	mat temp, vt1, vt2, vt3, vt4;
	temp = grad_slav;
	cout << "Successfully got Gradients from master of size: " << temp.n_rows << endl;
	Dl.col(7) = temp.rows(0, 166);
	Dl.col(9) = temp.rows(167, 333);
	Dl.col(11) = temp.rows(334, 500);
	vt1 = temp.rows(501, 503);
	vt2 = temp.rows(504, 506);
	vt3 = temp.rows(507, 509);
	vt4 = temp.rows(510, 512);
	//cout << "vt1: \n" << vt1 << endl;
	Dl(0,span(13,15)) = vt1.t();
	Dl(12,span(13,15)) = vt2.t();
	Dl(34,span(13,15)) = vt3.t();
	Dl(153,span(13,15)) = vt4.t();
	m_Synchronized = 0;
}*/

if(m_Synchronized1==1 || m_Synchronized2==1 || m_Synchronized3==1 || m_Synchronized4==1 || m_Synchronized5==1 || m_Synchronized6==1 || m_Synchronized7==1 || m_Synchronized8==1 || m_Synchronized9==1 || m_Synchronized10==1 || m_Synchronized11==1 || m_Synchronized12==1 || m_Synchronizedvr1==1 || m_Synchronizedvr2==1 || m_Synchronizedvr3==1 || m_Synchronizedvr4==1)
{
	t4 = boost::posix_time::second_clock::local_time();
	dt = t4 - t3;
	delay(l1) = (dt.minutes()*60) + dt.seconds();
	l1 = l1+1;
	delay.save("delay.csv",csv_ascii);
}

if (m_Synchronized1)
{
	Dl(19,7) = grad_slav1(0);
	Dl(22,7) = grad_slav1(0);
	cout << "Successfully got Gradients from slave 1 of size: " << grad_slav1.n_rows << endl;
	m_Synchronized1 = 0;	
}
else
{
	cout << "couldn't Load Gradients from slave 1" << endl;
	Dl(19,7) = Dl_old(19,7);
	Dl(22,7) = Dl_old(19,7);
}

if (m_Synchronized2)
{
	Dl(90,11) = grad_slav2(2);
	cout << "Successfully got Gradients from slave 2 of size: " << grad_slav2.n_rows << endl;
	m_Synchronized2 = 0;	
}
else
{
	cout << "couldn't Load Gradients from slave 2" << endl;
	Dl(90,11) = Dl_old(90,11);
}

if (m_Synchronized3)
{
	Dl(85,9) = grad_slav3(1);
	cout << "Successfully got Gradients from slave 3 of size: " << grad_slav3.n_rows << endl;
	m_Synchronized3 = 0;	
}
else
{
	cout << "couldn't Load Gradients from slave 3" << endl;
	Dl(85,9) = Dl_old(85,9);
}

if (m_Synchronized4)
{
	Dl(138,7) = grad_slav4(1);
	Dl(145,11) = grad_slav4(2);
	cout << "Successfully got Gradients from slave 4 of size: " << grad_slav4.n_rows << endl;
	m_Synchronized4 = 0;	
}
else
{
	cout << "couldn't Load Gradients from slave 4" << endl;
	Dl(138,7) = Dl_old(138,7);
	Dl(145,11) = Dl_old(145,11);
}

if (m_Synchronized5)
{
	Dl(49,7) = grad_slav5(0);
	Dl(49,9) = grad_slav5(1);
	Dl(49,11) = grad_slav5(2);
	Dl(128,7) = grad_slav5(0);
	Dl(128,9) = grad_slav5(1);
	Dl(128,11) = grad_slav5(2);
	cout << "Successfully got Gradients from slave 5 of size: " << grad_slav5.n_rows << endl;
	m_Synchronized5 = 0;	
}
else
{
	cout << "couldn't Load Gradients from slave 5" << endl;
	Dl(49,7) = Dl_old(49,7);
	Dl(49,9) = Dl_old(49,9);
	Dl(49,11) = Dl_old(49,11);
	Dl(128,7) = Dl_old(128,7);
	Dl(128,9) = Dl_old(128,9);
	Dl(128,11) = Dl_old(128,11);
}

if (m_Synchronized6)
{
	Dl(63,7) = grad_slav6(0);
	Dl(63,9) = grad_slav6(1);
	Dl(63,11) = grad_slav6(2);
	cout << "Successfully got Gradients from slave 6 of size: " << grad_slav6.n_rows << endl;
	m_Synchronized6 = 0;	
}
else
{
	cout << "couldn't Load Gradients from slave 6" << endl;
	Dl(63,7) = Dl_old(63,7);
	Dl(63,9) = Dl_old(63,9);
	Dl(63,11) = Dl_old(63,11);
}

if (m_Synchronized7)
{
	Dl(70,7) = grad_slav7(0);
	Dl(70,9) = grad_slav7(1);
	Dl(70,11) = grad_slav7(2);
	cout << "Successfully got Gradients from slave 7 of size: " << grad_slav7.n_rows << endl;
	m_Synchronized7 = 0;	
}
else
{
	cout << "couldn't Load Gradients from slave 7" << endl;
	Dl(70,7) = Dl_old(70,7);
	Dl(70,9) = Dl_old(70,9);
	Dl(70,11) = Dl_old(70,11);
}

if (m_Synchronized8)
{
	Dl(6,7) = grad_slav8(0);
	cout << "Successfully got Gradients from slave 8 of size: " << grad_slav8.n_rows << endl;
	m_Synchronized8 = 0;	
}
else
{
	cout << "couldn't Load Gradients from slave 8" << endl;
	Dl(6,7) = Dl_old(6,7);
}

if (m_Synchronized9)
{
	Dl(40,11) = grad_slav9(2);
	cout << "Successfully got Gradients from slave 9 of size: " << grad_slav9.n_rows << endl;
	m_Synchronized9 = 0;	
}
else
{
	cout << "couldn't Load Gradients from slave 9" << endl;
	Dl(40,11) = Dl_old(40,11);
}

if (m_Synchronized10)
{
	Dl(149,9) = grad_slav10(1);
	cout << "Successfully got Gradients from slave 10 of size: " << grad_slav10.n_rows << endl;
	m_Synchronized10 = 0;	
}
else
{
	cout << "couldn't Load Gradients from slave 10" << endl;
	Dl(149,9) = Dl_old(149,9);
}

if (m_Synchronized11)
{
	Dl(33,7) = grad_slav11(0);
	cout << "Successfully got Gradients from slave 11 of size: " << grad_slav11.n_rows << endl;
	m_Synchronized11 = 0;	
}
else
{
	cout << "couldn't Load Gradients from slave 11" << endl;
	Dl(33,7) = Dl_old(33,7);
}

if (m_Synchronized12)
{
	Dl(29,11) = grad_slav12(2);
	cout << "Successfully got Gradients from slave 12 of size: " << grad_slav12.n_rows << endl;
	m_Synchronized12 = 0;	
}
else
{
	cout << "couldn't Load Gradients from slave 12" << endl;
	Dl(29,11) = Dl_old(29,11);
}

if (m_Synchronizedvr1)
{
	Dl(0,13) = grad_slavvr1(0);
	Dl(0,14) = grad_slavvr1(1);
	Dl(0,15) = grad_slavvr1(2);
	cout << "Successfully got Gradients from slave vr1 of size: " << grad_slavvr1.n_rows << endl;
	m_Synchronizedvr1 = 0;	
}
else
{
	cout << "couldn't Load Gradients from slave vr1" << endl;
	Dl(0, span(13, 15)) = Dl_old(0, span(13, 15));
}

if (m_Synchronizedvr2)
{
	Dl(12,13) = grad_slavvr2(0);
	Dl(12,14) = grad_slavvr2(1);
	Dl(12,15) = grad_slavvr2(2);
	cout << "Successfully got Gradients from slave vr2 of size: " << grad_slavvr2.n_rows << endl;
	m_Synchronizedvr2 = 0;	
}
else
{
	cout << "couldn't Load Gradients from slave vr2" << endl;
	Dl(12, span(13, 15)) = Dl_old(12, span(13, 15));
}

if (m_Synchronizedvr3)
{
	Dl(34,13) = grad_slavvr3(0);
	Dl(34,14) = grad_slavvr3(1);
	Dl(34,15) = grad_slavvr3(2);
	cout << "Successfully got Gradients from slave vr3 of size: " << grad_slavvr3.n_rows << endl;
	m_Synchronizedvr3 = 0;	
}
else
{
	cout << "couldn't Load Gradients from slave vr3" << endl;
	Dl(34, span(13, 15)) = Dl_old(34, span(13, 15));
}

if (m_Synchronizedvr4)
{
	Dl(153,13) = grad_slavvr4(0);
	Dl(153,14) = grad_slavvr4(1);
	Dl(153,15) = grad_slavvr4(2);
	cout << "Successfully got Gradients from slave vr4 of size: " << grad_slavvr4.n_rows << endl;
	m_Synchronizedvr4 = 0;	
}
else
{
	cout << "couldn't Load Gradients from slave vr4" << endl;
	Dl(153, span(13, 15)) = Dl_old(153, span(13, 15));
}

VPQ dpf_re = DPF_return123(Dl, Z);
mat Vpolar = dpf_re.Vpolar;
mat PQb = dpf_re.PQb;
mat PQL = dpf_re.PQL;

int Lvp = Vpolar.n_rows;
int Wvp = Vpolar.n_cols;
int Lpq = PQb.n_rows;
int Wpq = PQb.n_cols;
double Pla_t = sum(PQL.col(0));
double Plb_t = sum(PQL.col(2));
double Plc_t = sum(PQL.col(4));

mat Pload_total, Ploss_orig_ph;
double Ploss_orig;
Pload_total<< Pla_t << Plb_t << Plc_t << endr;
Ploss_orig_ph << (PQb(0, 0) - Pla_t) << (PQb(0, 2) - Plb_t) << (PQb(0, 4) - Plc_t) << endr;
cout << "phase loss" << Ploss_orig_ph << endl;
Ploss_orig = accu(Ploss_orig_ph);

int Lnum_a = Y_return.Lnum_a;
int Lnum_b = Y_return.Lnum_b;
int Lnum_c = Y_return.Lnum_c;

Vabc V_abc = V_abc_list(Vpolar,Node_f,Lvp,Lnum_a,Lnum_b,Lnum_c);
mat V_a = V_abc.V_a;
mat V_b = V_abc.V_b;
mat V_c = V_abc.V_c;

mat Vmin_abc, Vmax_abc;
Vmin_abc << min(min(V_a)) << min(min(V_b)) << min(min(V_c))<<endr;
Vmax_abc << max(max(V_a)) << max(max(V_b)) << max(max(V_c))<<endr;
double Vmin_orig = min(min(Vmin_abc));
double Vmax_orig = max(max(Vmax_abc));

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

		/*BOOST_FOREACH(CPeerNode peer, m_peers | boost::adaptors::map_values)
        	{
      	    
	   	mat aVR1, aVR2, aVR3, aVR4;
		aVR1 = Dl(0, span(13, 15));
		aVR2 = Dl(12, span(13, 15));
		aVR3 = Dl(34, span(13, 15));
		aVR4 = Dl(153, span(13, 15));
		aVR1 = aVR1.t();
		aVR2 = aVR2.t();
		aVR3 = aVR3.t();
		aVR4 = aVR4.t();	
			
	    	mat S2;
		S2 = join_cols(join_cols(Dl.col(6), Dl.col(8)), Dl.col(10));
		S2 = join_cols(S2,Dl.col(7));
		S2 = join_cols(S2,Dl.col(9));
		S2 = join_cols(S2,Dl.col(11));
		S2 = join_cols(S2,aVR1);
		S2 = join_cols(S2,aVR2);
		S2 = join_cols(S2,aVR3);
		S2 = join_cols(S2,aVR4);
		std::cout << "S2: \n" << S2.n_rows << std::endl;
	        if (peer.GetUUID() == "arrogate:5001")
			{
				cout << "Sending to: " << peer.GetUUID() << " the Load Data!" << endl;
				ModuleMessage mg = Gradient(S2);
				peer.Send(mg);
			}	    
		}*/

//Time duration calculation
t2 = boost::posix_time::second_clock::local_time();

//t2(boost::posix_time::second_clock::local_time());
td = t2 - t1;
cout << "Time elapsed is: " << td.seconds() << endl;

time(l) = td.seconds();
l = l + 1;
time.save("time.csv",csv_ascii);

Dl_old = Dl;

if (td.minutes() > 1 || (td.minutes() == 1 && td.seconds() > 0) || (i==0 && td.seconds() > 53)	)
{
	min_v(i) = Vmin_orig;
	max_v(i) = Vmax_orig;
	total_loss(i) = Ploss_orig;

	min_v.save("min_volt.csv", csv_ascii);
	max_v.save("max_volt.csv", csv_ascii);
	total_loss.save("total_loss.csv", csv_ascii);
	
	//t5 = boost::posix_time::microsec_clock::local_time();
	//delt = t5 - t0;
	
	i = i + 1;
	t1 = t2;
	cout << "Moving on to next Load" << endl;
	Dl.col(6) = load_profile(i)*Dl_ld.col(6) - pv_profile(i)*Dl_pv.col(6);
	Dl.col(8) = load_profile(i)*Dl_ld.col(8) - pv_profile(i)*Dl_pv.col(8);
	Dl.col(10) = load_profile(i)*Dl_ld.col(10) - pv_profile(i)*Dl_pv.col(10);
	Dl.col(7) = load_profile(i)*(Dl_ld.col(7) - Dl_pv.col(7));
	Dl.col(9) = load_profile(i)*(Dl_ld.col(7) - Dl_pv.col(7));
	Dl.col(11) = load_profile(i)*(Dl_ld.col(7) - Dl_pv.col(7));
	for (int z = 0; z < 15; z++)
	{
		Dl(PV_node(i), 7) = Dl_old(PV_node(i), 7);
		Dl(PV_node(i), 9) = Dl_old(PV_node(i), 9);
		Dl(PV_node(i), 11) = Dl_old(PV_node(i), 11);
	}

	BOOST_FOREACH(CPeerNode peer, m_peers | boost::adaptors::map_values)
        	{
      	    
	   	mat aVR1, aVR2, aVR3, aVR4;
		aVR1 = Dl(0, span(13, 15));
		aVR2 = Dl(12, span(13, 15));
		aVR3 = Dl(34, span(13, 15));
		aVR4 = Dl(153, span(13, 15));
		aVR1 = aVR1.t();
		aVR2 = aVR2.t();
		aVR3 = aVR3.t();
		aVR4 = aVR4.t();	
			
	    	mat S2;
		S2 = join_cols(join_cols(Dl.col(6), Dl.col(8)), Dl.col(10));

		S2 = join_cols(S2,Dl.col(7));
		S2 = join_cols(S2,Dl.col(9));
		S2 = join_cols(S2,Dl.col(11));
		S2 = join_cols(S2,aVR1);
		S2 = join_cols(S2,aVR2);
		S2 = join_cols(S2,aVR3);
		S2 = join_cols(S2,aVR4);
		std::cout << "S2: \n" << S2.n_rows << std::endl;

	        if (peer.GetUUID() == "arrogate:5001")
			{
				cout << "Sending to: " << peer.GetUUID() << " the Load Data!" << endl;
				ModuleMessage mg = Gradient(S2);
				peer.Send(mg);
			}	    

		}
}


}// end of vvc_main()


}//namespace vvc
}// namespace broker
}// namespace freedm
//}
