import React, { Component } from "react";
//import { Spring } from "react-spring/renderprops";
import UserTranscripts from "./ChildComponents1level/UserTranscripts";
import AllTranscripts from "./ChildComponents1level/AllTranscripts";
import { MDexonalerts, MDlastexon, MDgeneCard, MDexoncodes, MDaddNightingaleReadyJSON } from "./ModuleHF";
export class HolderFront extends Component {
  constructor(props) {
    super(props);
    this.state = {
      scrollable: false,
      showTrans: true,
      showUserTrans: false
    };
  }
  toggleTranscript = e =>
    this.setState({
      showTrans: true,
      showUserTrans: false
    });
  toggleUserTrans = e =>
    this.setState({
      showTrans: false,
      showUserTrans: true
    });
  render() {
    // let gene = require("../oldSet/sample.json");
    
    let gene= this.props.gene;
    // console.log('JSON',gene)
    let exonHash = {};
    let exonalertscodes = MDexonalerts(gene.exonsGenes);
    // console.log("before this1");
    gene.exonsGenes.map(ex => (exonHash[ex.id] = ex));
    // console.log(exonHash)
    // console.log("there it is0")
    gene=MDaddNightingaleReadyJSON(gene,exonalertscodes,exonHash)
    // console.log("there it is1")

    //console.log("in holderfront");
    //console.log("gene", gene);
    // console.log("transcriptsall", this.state.showTrans);
    // console.log("transcriptUser", this.state.showUserTrans);


    /////////////////////////////////////////////////////////////////
    /// secondary structure addition will be default              ///
    /// decide the shape of the disoredr region                   ///
    /// domains will be 1 interpro track and exonw will be another///
    /// also have to deign the datatable, that will enjoy the exon///
    /// alert codes                                               ///


    //iterate gene JSON, for everyofits exonsGenes keys, get its internal sql id and store exon
    let showtransclass = this.state.showTrans ? "nav-link active text-primary" : "nav-link text-secondary";
    let showusertransclass = this.state.showUserTrans
      ? "nav-link active text-primary"
      : "nav-link text-secondary";
    console.log("before this2");
    let exoncod = MDexoncodes(gene.exonsGenes);
    console.log("before this3s");
    return (
      <div className="container  pb-5 mb-5">
        {/* style={{borderStyle: "dotted",  borderColor: "blue"}} */}
        <div className="pb-2 m-auto">
            {MDgeneCard(gene)}
        </div>

          <div
            className="row justify-content-center mt-2 pt-2"
            style={{ paddingTop: "0rem", borderCollapse: "collapse" }}>
              <div className="col md-10">
                  <ul className="nav nav-tabs ">
                    <li className="nav-item">
                      <a
                        className={showtransclass}
                        href="#"
                        onClick={this.toggleTranscript}
                        style={{ fontSize: 20, fontWeight: 500 }}
                        
                      >
                        Show transcripts
                      </a>
                    </li>
                    <li className="nav-item">
                      <a
                        onClick={this.toggleUserTrans}
                        className={showusertransclass}
                        href="#"
                        style={{ fontSize: 20, fontWeight: 500 }}
                      >
                        Enter your protein sequence
                      </a>
                    </li>
                  </ul>
            </div>
          </div>
        <div
          className="row justify-content-center"
          style={{ paddingTop: "0.5rem" }}
        >
          <div
            className="col"
            style={{ display: "inline-block", alignContent: "center" }}
          >
            {this.state.showTrans ? (
              <AllTranscripts
                transSet={gene.trans}
                domset={gene.domainCod}
                exOn={exonHash}
                finalExon={MDlastexon(gene.exonsGenes)}
                exonCodes={exoncod}
                exonAlerts={exonalertscodes}
              />
            ) : null}
            {this.state.showUserTrans ? (
              <UserTranscripts
                transet={gene.trans}
                exOn={gene.exonsGenes}
                exonCodes={exonHash}
                exonAlerts={exonalertscodes}
              />
            ) : null}
          </div>
        </div>
      </div>
    )
  }
  // else take the object in form of ssor domain and render below
}
export default HolderFront;
