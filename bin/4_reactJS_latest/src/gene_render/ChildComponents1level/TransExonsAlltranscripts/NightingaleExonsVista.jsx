import React, { Component, Fragment, useState, useEffect, useCallback } from "react";
import ProtvistaManager from "protvista-manager";
import ProtvistaTrack from "protvista-track";
import ProtvistaNavigation from "protvista-navigation";
import ProtvistaSequence from "protvista-sequence";
import ProtvistaInterproTrack from "protvista-interpro-track";
import ProtvistaColouredSequence from "protvista-coloured-sequence";
import loadWebComponent from "./utils/load-web-component";
//import sequence from "../mocks/sequence.json";
// import { dataIPR, signatures } from "../mocks/interpro";
// import secondaryStructureData from "../mocks/interpro-secondary-structure.json";
import ProtvistaSaver from "protvista-saver";
import ProtvistaOverlay from "protvista-overlay";
import ProtvistaZoomTool from "protvista-zoom-tool";
import ProtvistaDatatable from "protvista-datatable";
//import { dataTable } from "../mocks/datatable";
import Collapsible from 'react-collapsible';
import Accordion from 'react-bootstrap/Accordion'
/*
Three tracks, exons, domains, secondary structres and disorder, exonTrack
*/
class NightingaleExonsVista extends Component {
  constructor(props){
    super(props);
    //console.log ('props',this.props)
    this.state = {
      sequence: this.props.sequence,
      domains: this.props.nightDom,
      ss: this.props.nightSs,
      exons: this.props.nighthTable,
  }}
  componentDidMount() {
     console.log("cdm1",this.props.transwhole)
     console.log("cdm2",this.props.transwhole.nightDom)
     //this["interpro-track"].innerHTML = this.props.transwhole.nightDom;
     let seqTrack="#sequence-track" + this.props.variable
     console.log(seqTrack)
     document.querySelector("#sequence-track"+this.props.variable).data = this.props.transwhole.sequence;
     document.querySelector("#domain-track"+this.props.variable).data = this.props.transwhole.nightDom;
     document.querySelector("#exon-track"+this.props.variable).data = this.props.transwhole.nightExons;
     document.querySelector("#secondaryStructure-track"+this.props.variable).data = this.props.transwhole.nightSs;
     document.querySelector("#disorder-track"+this.props.variable).data = this.props.transwhole.nightDis;     

    // document.querySelector("#saver2"+this.props.variable).backgroundColor = "#ddddee";
    }

    render () {
    loadWebComponent("protvista-manager", ProtvistaManager);
    loadWebComponent("protvista-track", ProtvistaTrack);
    loadWebComponent("protvista-navigation", ProtvistaNavigation);
    loadWebComponent("protvista-sequence", ProtvistaSequence);
    loadWebComponent("protvista-interpro-track", ProtvistaInterproTrack);
    loadWebComponent("protvista-saver", ProtvistaSaver);
    loadWebComponent("protvista-overlay", ProtvistaOverlay);
    loadWebComponent("protvista-zoom-tool", ProtvistaZoomTool);
    loadWebComponent("protvista-datatable", ProtvistaDatatable);
    console.log("render",this.props.nightDom)
    console.log("render2",this.props)
    console.log("STATE",this.state)
    console.log("TABLE", this.props.transwhole.nightExons)
    const lengthSeq=this.props.transwhole.sequence.length
    return (
      <div className="row px-2 mx-2 my-2 py-1" style = {{borderStyle:"double",  borderWidth: "0 0 5px 0" /* top right bottom left */ }}>
          <Fragment>
            {/* <protvista-saver
              element-id="just-tracks" id={"saver2"+this.props.variable} file-name="tracks" file-format="jpeg">
              <button>Download Tracks</button>
            </protvista-saver> */}
            <protvista-manager
              displaystart= "1" displayend="100" length ={lengthSeq}>
              <div className="row border py-0 my-0 px-2 mx-2 align-items-center">
                <div className="col-2 px-2 mx-4"style={{float:"left",
                  // color: "white",
                  // borderRadius: "4px",
                  // margin: "2px",
                  // verticalAlign:"middle"
                 }}>
                    <span className="badge align-center mx-2 my-3 h5"
                      style = {{color: "#00303f", background: "#e9e9e9", border:"1px solid #00303f", borderLeft: "5px solid #00303f"}}
                    >{this.props.transwhole.tId}
                    </span>
                    {/* <span>{this.props.transwhole.tId}</span> */}
                </div>
                <div className="col-5"></div>
                <div className="col-4 mx-auto px-auto align-self-center" align="right">
                  <protvista-zoom-tool
                    length={lengthSeq}
                    displaystart="1"
                    displayend="100"
                    style={{
                      fontSize: "medium",
                      "--button-font-size" : "12px",
                      "--button-padding": "0px",
                      "--button-margin": "0px",
                      "--button-background": "#a7d2cb",
                      "--button-text-color": "#00303f",
                      "--button-background-focus": "#00a6d5",
                      "--button-border-radius": "5px"}}>
                      <span className="m-0 p-0" slot="zoom-in">+</span>
                      <span className="m-0 p-0" slot="zoom-in-seq">Zoom to Sequence</span>
                  </protvista-zoom-tool>
                </div>
              </div>
              
              <div className="row p-3 my-1 px-2 mx-2">
              <protvista-navigation length ={lengthSeq}  /> 
                  <div id={"just-tracks"}>
                      <protvista-sequence
                        length={lengthSeq} displaystart="1" displayend="100"
                        id={"sequence-track"+ this.props.variable}
                        highlight-event="onmouseover"
                        use-ctrl-to-zoom />
                      <protvista-track
                          id={"exon-track"+ this.props.variable}
                          length={lengthSeq}
                          // layout="non-overlapping"
                          highlight-event="onmouseover"
                          expanded
                          use-ctrl-to-zoom
                        />
                      <protvista-track
                          id={"secondaryStructure-track" + this.props.variable}
                          length={lengthSeq}
                          layout="non-overlapping"
                          use-ctrl-to-zoom
                        />
                        <protvista-track
                          id={"disorder-track" + this.props.variable}
                          length={lengthSeq}
                          layout="non-overlapping"
                          use-ctrl-to-zoom
                        />
                        <protvista-interpro-track
                          id={"domain-track" + this.props.variable}
                          length={lengthSeq}
                          highlight-event="onmouseover"
                          expanded
                          use-ctrl-to-zoom
                        />
                    </div>
                </div>
                <div className="row mx-0 my-0 py-0">
                <Accordion flush>
                  <Accordion.Item eventKey="0">
                    <Accordion.Header>
                      <div className="badge px-3 py-2" style = {{color: "#00303f", fontSize:"15px", fontWeight:"normal", background: "#a7d2cb"}}>
                      Show/Hide Details
                      </div>
                    </Accordion.Header>
                  <Accordion.Body>
                    {/* <h5 style={{paddingBottom:"20px"}}>DataTable</h5> */}
                  <protvista-datatable displaystart="1" displayend="100">
                    <table>
                      <thead>
                        <tr>
                          <th data-filter="ft_key">Feature key </th>
                          <th>Description </th>
                          <th>Positions </th>
                        </tr>
                      </thead>
                      <tbody>
                        {/* {console.log("dataNew",dataTable)} */}
                        {this.props.transwhole.nighthTable?.map((row, i) => (
                          //i is row index, starting from 0
                          <Fragment key={i}>
                            {console.log("da*",row,i)}
                            <tr
                              data-id={`${row.start}-${row.end}`}
                              data-start={row.start}
                              data-end={row.end}
                            >
                              {/* inan loop and row content is decalred */}
                              <td data-filter="ft_key" data-filter-value={row.type}>
                                {row.type}
                              
                              </td>
                              <td>{row.description}</td>
                              <td>
                                {row.start} - {row.end}
                              </td>
                            </tr>
                            <tr data-group-for={`${row.start}-${row.end}`}>
                              <td>{row.tooltip}</td>
                            </tr>
                          </Fragment>
                        ))}
                      </tbody>
                    </table>
                    
                  </protvista-datatable>
                  </Accordion.Body>
                  </Accordion.Item>
                  </Accordion>
                  </div>
                  
                
        </protvista-manager>
      </Fragment>
      </div>
    )
  }
  }

export default NightingaleExonsVista;
