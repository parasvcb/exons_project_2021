import React, { Component } from "react";
import ChooseTranscripts from "./TransExonsAlltranscripts/ChooseTranscripts";
import RenderExonsFull from "./TransExonsAlltranscripts/RenderExonsFull_vista";
import RenderExonsAln from "./TransExonsAlltranscripts/RenderExonsAln";
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
    return (
      <>
        <ChooseTranscripts
          trans={transcript}
          handleUpdate={this.handleClick}
          transadded={this.state.arraylis}
        />

        {this.state.arraylis.length > 0 ? (
          <div className="container-fluid align-items-center">
            <ul className="nav nav-tabs ">
              <li className="nav-item">
                <a href="#" onClick={this.toggle_exonfull}>
                  <span
                    className={this.state.exonfull ? "badge badge-dark" : "badge badge-light"}
                    style={{fontSize: this.state.exonfull ? 20 : 15, fontWeight: 500 }}
                    //color: this.state.exonfull ? "white" : "black"
                  >
                    Exon details
                  </span>
                </a>
              </li>
              <li className="nav-item">
                <a onClick={this.toggle_exonNightingale} href="#">
                  <span
                    className={this.state.exonaln ? "badge badge-dark" : "badge badge-light" }
                    style={{ fontSize: this.state.exonNightingale ? 20 : 15, fontWeight: 500 }}
                      //</a>color: this.state.exonNightingale ? "white" : "black"
                  >
                    Exon Nightingale
                  </span>
                </a>
              </li>
              <li className="nav-item">
                <a onClick={this.toggle_exonaln} href="#">
                  <span
                    className={this.state.exonaln ? "badge badge-dark" : "badge badge-light"}
                    style={{ fontSize: this.state.exonaln ? 20 : 15, fontWeight: 500 }}
                      //color: this.state.exonaln ? "white" : "black"
                  >
                    Exon alignment
                  </span>
                </a>
              </li>
            </ul>
            {this.state.exonfull ? (
              <>
                <button
                  className="m-2 btn btn-warning btn-md"
                  onClick={this.handleExonInf}
                >
                  {this.state.exoninf ? "Hide exons nomenclature" : "Show exons nomenclature"}
                </button>
                <div className="row">
                  <div className="col-sm-9 border">
                    {this.state.arraylis.map(tr1 => (
                      <div key={tr1.tId}>
                        <RenderExonsFull
                          transwhole={tr1} domset={this.props.domset} handleUpdateInfRend={this.showInf} addAlert={this.addAlert} removeAlert={this.removeAlert} exOn={this.props.exOn} exonAlerts={this.props.exonAlerts}
                        />
                      </div>
                    ))}
                  </div>
                  <div className="col-sm-3 border sticky-top">
                    <div className="d-flex flex-column sticky-top">
                      <div
                        className="p-2 bg-info"
                        style={{
                          Height: "350px"
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
            ) : this.state.exonNightingale ? (<h1>test</h1>) :
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

export default AllTranscripts;

/*
This componnent will send the data of transcripts to the choose_transcripts comonent
that component will add data to state of this component
and that state from this will again be sent to the exon comp to be renedered
therw ill be two exons copmonents btw

*/
