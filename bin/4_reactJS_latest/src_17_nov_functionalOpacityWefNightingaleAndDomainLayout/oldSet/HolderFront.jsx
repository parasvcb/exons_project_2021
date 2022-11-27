import React, { Component } from "react";
import UserTranscripts from "./ChildComponents1level/UserTranscripts";
import AllTranscripts from "./ChildComponents1level/AllTranscripts";
import { memberExpression } from "@babel/types";

const gene = require("./103.json");

export class HolderFront extends Component {
  constructor(props) {
    super(props);
    this.state = {
      gene,
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
    console.log("gene", this.state.gene);
    console.log("transcriptsall", this.state.showTrans);
    console.log("transcriptUser", this.state.showUserTrans);
    let exonalertscodes = this.exonalerts(this.state.gene.exonsGenes);
    let exonHash = {};
    this.state.gene.exonsGenes.map(ex => (exonHash[ex.id] = ex));
    let showtransclass = this.state.showTrans ? "nav-link active" : "nav-link";
    let showusertransclass = this.state.showUserTrans
      ? "nav-link active"
      : "nav-link";
    let exoncod = this.exoncodes(this.state.gene.exonsGenes);
    return (
      <div className="container-fluid">
        {this.geneCard()}
        <div
          className="row justify-content-center"
          style={{ paddingTop: "0rem", borderCollapse: "collapse" }}
        >
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
                transSet={this.state.gene.trans}
                domset={this.state.gene.domainCod}
                exOn={exonHash}
                finalExon={this.lastexon(this.state.gene.exonsGenes)}
                exonCodes={exoncod}
                exonAlerts={exonalertscodes}
              />
            ) : null}
            {this.state.showUserTrans ? (
              <UserTranscripts
                transet={this.state.gene.trans}
                exOn={this.state.gene.exonsGenes}
                exonCodes={exonHash}
              />
            ) : null}
          </div>
        </div>
      </div>
    );
  }
  lastexon(exons) {
    exons.sort(function(a, b) {
      return a.rawst - b.rawst;
    });
    const lastele = exons[exons.length - 1];
    if (lastele != "R") {
      return parseInt(lastele.exId.split(".")[3]);
    } else {
      return parseInt(lastele.exId.split(":")[2].split(".")[3]);
    }
  }
  exonalerts(exons) {
    let ob = {};
    exons.forEach(element => {
      let messageComp1 = "";
      if (element.exId[0] !== "R") {
        const id = element.exId.split(".");
        let message1h = id[0] + "." + id[1];
        let message1 = "";
        switch (id[0]) {
          case "U":
            message1 +=
              "This exon is always part of UTR in all the isoforms(ISFs) and always non-coding ";
            break;
          case "T":
            message1 +=
              "This exon is always part of CDS in all the isoforms(ISFs) ";
            break;
          case "M":
            message1 +=
              "This exon is always part of CDS in all the isoforms(ISFs) but has only 1nt contribution and hence no Amino acid ";
            break;
          case "D":
            message1 +=
              "This exon is either part of CDS or UTR in different isoforms(ISFs) of this gene ";
            break;
          default:
            message1 += " ";
        }
        if (
          ["D-1", "D-2", "D0", "T-1", "T-2", "T0"].indexOf(id[0] + id[1]) > -1
        ) {
          switch (id[0] + id[1]) {
            case "D-1":
              message1 +=
                "and it is non-coding due to PTC in upstream region for this ISF;";
              break;
            case "D-2":
              message1 += "and it is non-coding in this ISF;";
              break;
            case "D0":
              message1 += "and it has only 1nt share to CDS in this ISF;";
              break;
            case "T-1":
              message1 +=
                "but is not coding for aa in this ISF due PTC in upstream region;";
              break;
            case "T0":
              message1 += "but it has only 1nt share to CDS for this ISF;";
              break;
          }
        } else {
          switch (id[1]) {
            case "-2":
              message1 += "";
              break;
            case "-1":
              message1 += "";
              break;
            case "0":
              message1 += "";
              break;
            case "1":
              message1 += " and is coding amino acids (Reference) in this ISF,";
              break;
            default:
              message1 +=
                " and is coding amino acids different than that of reference exon (Variant:" +
                id[1] +
                ") in this ISF,";
          }
        }

        let message2h = id[2] + "." + id[3];
        let message2 =
          " further this exon is " + id[3] + " in gene sequence rank and ";
        switch (id[2]) {
          case "G":
            message2 += " is constitutively present in all the ISFs,";
            break;
          case "A":
            message2 +=
              " is alternatively present in ISFs with WEF:" + element.wef + ",";
            break;
          case "F":
            message2 +=
              " is present in all the ISFs with certain alternative splice site variations,";
            break;
          default:
            message2 += "";
        }

        let message3h = id[4] + "." + id[5];
        let message3 = "";
        switch (id[4]) {
          case "0":
            message3 +=
              " no alternative splice is choosen for this exon in this ISF";
            break;
          case "n":
            message3 +=
              " alternative 5' splice site has been choosen for this exon different from reference exon ";
            break;
          case "c":
            message3 +=
              " alternative 3' splice site has been choosen for this exon different from reference exon ";
            break;
          case "b":
            message3 +=
              " alternative 5' splice site and alternative 3' splice site has been choosen which are different from reference exon ";
            break;
          default:
            message3 += "";
        }
        if (["n", "c", "b"].indexOf(id[4]) > -1) {
          message3 += element.parent;
          message3 += " and that is " + id[5] + "such occurence";
        }
        messageComp1 = (
          <div>
            <p style={{ display: "inline", fontWeight: "bold" }}>
              <span style={{ color: "crimson" }}>{id[0] + "." + id[1]}.</span>
              <span style={{ color: "DarkGreen" }}>{id[2] + "." + id[3]}.</span>
              <span style={{ color: "DarkOrange" }}>
                {id[4] + "." + id[5]}.
              </span>
            </p>

            <p>
              <span style={{ color: "crimson", fontWeight: "bold" }}>
                ({message1h}){" "}
              </span>
              <span>{message1}</span>
              <span>{message2}</span>
              <span style={{ color: "DarkGreen", fontWeight: "bold" }}>
                ({message2h}){" "}
              </span>
              <span>{message3}</span>
              <span style={{ color: "DarkOrange", fontWeight: "bold" }}>
                ({message3h}){" "}
              </span>
            </p>
          </div>
        );
      } else {
        let id = element.exId.split(":");
        let message1h = id[0] + ":" + id[1];
        let message2h = id[2];
        let message3h = id[3];
        let message4h = id[4];
        let message1 = "This is a intron retention case which";
        switch (id[1]) {
          case "-2":
            message1 += " is non-coding in this transcript";
            break;
          case "1":
            message1 += " is coding in this transcript ";
            break;
          case "0":
            message1 += " has only 1nt contribution and hence no amino acid ";
            break;
          case "-1":
            message1 +=
              " is not coding for aa in this transcript due PTC in upstream region ";
            break;
          default:
            message1 += " ";
        }
        let message2 = " and spans from exon ";
        let message3 = " and is the ";
        let message4 = " and ending at exon ";
        messageComp1 = (
          <div>
            <h5>
              <p style={{ display: "inline", fontWeight: "bold" }}>
                <span style={{ color: "crimson" }}>{id[0] + ":" + id[1]}:</span>
                <span style={{ color: "DarkGreen" }}>{id[2]}:</span>
                <span style={{ color: "DarkOrange" }}>{id[3]}:</span>
                <span style={{ color: "SteelBlue" }}>{id[4]}</span>
              </p>
            </h5>
            <p>
              <span style={{ color: "crimson", fontWeight: "bold" }}>
                {message1h}{" "}
              </span>
              <span>{message1}</span>
              <span>{message2}</span>
              <span style={{ color: "DarkGreen", fontWeight: "bold" }}>
                {id[2]}
              </span>
              <span>{message3}</span>
              <span style={{ color: "DarkOrange", fontWeight: "bold" }}>
                {id[3]}
              </span>
              <span> instance starting from this exon</span>
              <span>{message4}</span>
              <span style={{ color: "SteelBlue", fontWeight: "bold" }}>
                {id[4]}
              </span>
            </p>{" "}
          </div>
        );
      }
      ob[element.exId] = messageComp1;
    });
    return ob;
  }

  exoncodes(exons) {
    let finalob = {};
    exons.forEach(element => {
      if (element.exId[0] !== "R" && element.exId.split(".")[4] === "0") {
        finalob[parseInt(element.exId.split(".")[3])] = element.exId[0];
      }
    });
    console.log("EXONCODES", finalob);
    return finalob;
  }
  geneCard() {
    const { gene } = this.state;
    let cod = 0;
    let non = 0;
    gene.exonsGenes.forEach(element => {
      if (element.aaseq.length > 0) {
        cod = cod + 1;
      } else {
        non = non + 1;
      }
    });
    return (
      <div className="d-flex">
        <div className="p-2 w-75 bg-light">
          <div className="card" style={cardStyle1}>
            <a
              target="_blank"
              href={"https://www.ncbi.nlm.nih.gov/gene/" + gene.entrezid}
            >
              <h4 className="card-title m-0 p-0 text-white">
                Gene: {gene.name}
              </h4>
            </a>
            <div className="card-body p-2 mb-2 text-white">
              <a
                target="_blank"
                href={"https://www.ncbi.nlm.nih.gov/gene/" + gene.entrezid}
              >
                <p className="m-0 p-0 text-white">NCBI ID: {gene.entrezid}</p>
              </a>
              <a
                target="_blank"
                href={
                  "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=" +
                  gene.txid
                }
              >
                <p className="m-0 pb-2 text-white">Organism: {gene.organism}</p>
              </a>
            </div>
          </div>
        </div>
        <div className="p-2 w-25 bg-light  ">
          <div className="card" style={cardStyle2}>
            <h4 className="card-title m-0 p-0">Constituents</h4>
            <div className="card-body p-1 m-0">
              <p className="m-0 p-0">
                Transcripts: {this.state.gene.trans.length}
              </p>
              <p className="m-0 p-0">Coding Exons: {cod}</p>
              <p className="m-0 p-0">NonCoding Exons: {non}</p>
            </div>
          </div>
        </div>
      </div>
    );
  }

  // else take the object in form of ss or domain and render below
}
const cardStyle1 = {
  background: "steelblue",

  paddingBottom: "10",
  color: "white"
};
const cardStyle2 = {
  background: "dodgerblue",
  color: "white",

  paddingBottom: "10"
};
export default HolderFront;
