import React, { Component } from "react";
import ChooseTranscripts from "./TransExonsAlltranscripts/ChooseTranscripts";
import RenderExonsFull from "./TransExonsAlltranscripts/RenderExonsFull_vista";
import RenderExonsAln from "./TransExonsAlltranscripts/RenderExonsAln";
import NightingaleExonsVista from "./TransExonsAlltranscripts/NightingaleExonsVista";
import {showInfFunc, showalert} from "./ModuleAT";
export class AllTranscripts extends Component {
  constructor(props) {
    super(props);
    this.state = {
      arraylis: [],
      exonfull: true,
      exonaln: false,
      exonNightingale:false,
      scrollable: false,
      alert: false,
      exoninf: true
    };
    this.handleClick = this.handleClick.bind(this);
    this.showInf = this.showInf.bind(this);
    this.addAlert = this.addAlert.bind(this);
    this.removeAlert = this.removeAlert.bind(this);
    //added state properties to them
  }
  toggle_exonfull = e =>
    this.setState({
      exonfull: true,
      exonNightingale:false,
      exonaln: false
    });
    toggle_exonNightingale = e =>
    this.setState({
      exonfull: false,
      exonNightingale:true,
      exonaln: false
    });
  toggle_exonaln = e =>
    this.setState({
      exonfull: false,
      exonNightingale:false,
      exonaln: true
    });
  addAlert(exid) {
    this.setState({
      alert: exid
    });
  }
  removeAlert() {
    this.setState({
      alert: false
    });
  }
  handleExonInf = e =>
    this.setState({
      exoninf: !this.state.exoninf
    });
  handleClick(tidinp) {
    let arr = [...this.state.arraylis];
    //console.log("indexof", arr.indexOf(tidinp), tidinp, arr);
    let ind = arr.indexOf(tidinp);
    //console.log("ind", ind);
    // let ind = "abc";

    if (ind !== -1) {
      arr.splice(ind, 1);
      this.setState({
        arraylis: arr
      });
    } else {
      this.setState({
        arraylis: [...this.state.arraylis, tidinp]
      });
    }
  }
  showInf(value) {
    this.setState({
      scrollable: value
    });
  }

  render() {
    const transcript = this.props.transSet;
    const finallist=this.state.arraylis[this.state.arraylis.length-1]
    let finalelem=finallist
    console.log("state",this.state)
    // console.log('&',typeof(finallist))    
    // if (this.state.arraylis.length>1){
    //    console.log(this.state.arraylis.length)
    //  console.log("finalelem2",finalelem.tId)
    //  console.log("finalelem3",finalelem)
    
    // }
    return (
      <>
        <ChooseTranscripts
          trans={transcript}
          handleUpdate={this.handleClick}
          transadded={this.state.arraylis}
        />

        {this.state.arraylis.length > 0 ? (
          <div className="container-fluid align-items-center mt-5">
            <ul className="nav nav-tabs ">
              <li className="nav-item">
                <a href="#" className="text-decoration-none" onClick={this.toggle_exonfull}>
                  <span
                    className={this.state.exonfull ? "nav-link active" : "nav-link"}
                    style={this.state.exonfull ? active3CompExondetailOn : active3CompExondetailOff}
                    // style={{fontSize: this.state.exonfull ? 20 : 15, fontWeight: 500}}
                    //color: this.state.exonfull ? "white" : "black"
                  >
                    Exon Details
                  </span>
                </a>
              </li>
              <li className="nav-item">
                <a onClick={this.toggle_exonNightingale} href="#" className="text-decoration-none">
                  <span
                    className={this.state.exonNightingale ? "nav-link active" : "nav-link"}
                    style={this.state.exonNightingale ? active3CompExondetailOn : active3CompExondetailOff}
                   
                      //</a>color: this.state.exonNightingale ? "white" : "black"
                  >
                    Exon Nightingale
                  </span>
                </a>
              </li>
              <li className="nav-item">
                <a onClick={this.toggle_exonaln} href="#" className="text-decoration-none">
                  <span
                    className={this.state.exonaln ? "nav-link active" : "nav-link"}
                    style={this.state.exonaln ? active3CompExondetailOn : active3CompExondetailOff}
                      //color: this.state.exonaln ? "white" : "black"
                  >
                    Exon Alignment
                  </span>
                </a>
              </li>
            </ul>
            {this.state.exonfull ? (
              <>
                <div className="clearfix">
                    <button
                      className="btn"
                      style={this.state.exoninf ? cardExonInfOn : cardExonInfOff}
                      onClick={this.handleExonInf}
                    >
                      {console.log('parent',this.state.arraylis)}
                      {this.state.exoninf ? "Hide exons nomenclature" : "Show exons nomenclature"}
                    </button>

                </div>
                <div className="row">
                  <div className="col-sm-9">
                    {this.state.arraylis.map(tr1 => (
                      <div key={tr1.tId}>
                        <RenderExonsFull
                          transwhole={tr1} domset={this.props.domset} handleUpdateInfRend={this.showInf} addAlert={this.addAlert} removeAlert={this.removeAlert} exOn={this.props.exOn} exonAlerts={this.props.exonAlerts}
                        />
                      </div>
                    ))}
                  </div>
                  <div className="col-sm-3 sticky-top">
                    <div className="d-flex flex-column sticky-top">
                      <div
                        className="p-2"
                        //  bg-info
                        style={{
                          Height: "350px",
                          background: "#284b63"
                        }}
                      >
                        {showalert(this.state)}
                      </div>
                      <div className="p-2 bg-light">
                        {showInfFunc(this.state.scrollable, this.props)}
                      </div>
                    </div>
                  </div>
                </div>
              </>
            ) : this.state.exonNightingale ? (
              <div className="container h-auto">
                {this.state.arraylis.map(tr1 => (
                <div key={tr1.tId} className="container" style={{margin: "0 0 5px 0"}}>
                  <NightingaleExonsVista transwhole={tr1} variable={tr1.tId.split('.')[0]} />
                </div>
                  ))} 
              </div>
            ) :
              this.state.exonaln ? (
              <RenderExonsAln
                arraylis={this.state.arraylis} total={this.props.finalExon}exOn={this.props.exOn} exonCodes={this.props.exonCodes}
              />
            ):''} 
            </div> ): null}
      </>
    );
  }

  
}

/*
This program will be modified will fork into two components
one will hold the card and modify the state
the other will render the exons on basis of state of that values
create a list of states for each variable and then creat a function that modifies that state
once clicked in another compinenet this methdo shold be passed and caorresodong state shuld be canged
with this functin whch s apssed as prop

*/

const c1Style = {
  background: "steelblue",
  color: "white",
  padding: " 1.5rem"
};

const active3CompExondetailOn = {
  fontSize: 20, 
  fontWeight: "normal",
  color: "#dc3545",
  color: "#b37d4e"
}

const active3CompExondetailOff = {
  fontSize: 15, 
  fontWeight: "normal",
  color: "#00303f"
}

const cardExonInfOn = {
  float:"right",
  fontFamily: "Ubuntu Mono", 
  marginRight: "30px", 
  background:"none", 
  color: "#b37d4e",
  fontWeight: "bold",
  border: "double 1px black"
}

const cardExonInfOff = {
  float:"right", 
  fontFamily: "Ubuntu Mono", 
  marginRight: "30px", 
  background:"none", 
  color: "",
  border: "double 1px black"
}
// const active3CompExondetailOn = {
//   fontSize: this.state.exonfull ? 20 : 15, fontWeight: 500}}
// }
export default AllTranscripts;

/*
This componnent will send the data of transcripts to the choose_transcripts comonent
that component will add data to state of this component
and that state from this will again be sent to the exon comp to be renedered
therw ill be two exons copmonents btw

*/
