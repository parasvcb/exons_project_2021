import React, { Component } from "react";
import About from "./about";
import Tutorial from "./tutorial";
//import AppLand from "../landing_api/Landing_page";
import Nomenclature from "./nomenclature";
//import Badge from 'react-bootstrap/Badge';
import HolderFront from '../gene_render/HolderFront';

class Navbar extends Component {
  constructor(props) {
    super(props);
    this.state = {
      home: true,
      about: false,
      tutorial: false,
      nomenclature: false
    };
  }
  toggle_home = e =>
    this.setState({
      about: false,
      tutorial: false,
      nomenclature: false,
      home: true
    });

  toggle_name = e =>
    this.setState({
      home: false,
      about: false,
      tutorial: false,
      nomenclature: true
    });

  toggle_about = e =>
    this.setState({
      home: false,
      tutorial: false,
      nomenclature: false,
      about: true
    });

  toggle_tutorial = e =>
    this.setState({
      home: false,
      about: false,
      nomenclature: false,
      tutorial: true
    });

  render() {
    return (
      <>
      {/*console.log(this.state)*/}
        <nav className="nav" style={headerStyle}>
          <a
            className="nav-link active"
            //href="http://14.139.227.206/nextrap"
            href="http://localhost:3000"

            onClick={this.toggle_home}
            style={{ color: "white" }}
          >
            <span>Home</span>
          </a>
          <a
            className="nav-link"
            href="#"
            onClick={this.toggle_about}
            style={{ color: "white" }}
          >
            About
          </a>
          <a
            className="nav-link"
            href="#"
            onClick={this.toggle_name}
            style={{ color: "white" }}
          >
            Exon-nomenclature
          </a>
          <a
            className="nav-link"
            href="#"
            onClick={this.toggle_tutorial}
            style={{ color: "white" }}
          >
            Tutorial
          </a>
        </nav>
        <>
        {this.state.about ? <About/> : ' '}
        {this.state.nomenclature ? <Nomenclature/> : ' '}
        {/* {this.state.home ? <AppLand/> : ' '} */}
        {this.state.home ? <HolderFront/> : ' '}
        {this.state.tutorial ? <Tutorial/> : ' '}
        {console.log('InThere')}
        </>
      </> 
    );
  }
}
const headerStyle = {
  background: "Seagreen",
  fontWeight: 500,
  color: "Honeydew",
  padding: "0rem",
  fontFamily: "font-family: Arial, Helvetica, sans-serif",
  marginTop: "0px",
  marginBottom: "1.5rem"
};
const boldStyle = {
  color: "white",
  fontWeight: "thick"
};

export default Navbar;
