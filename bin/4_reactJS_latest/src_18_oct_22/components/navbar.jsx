import React, { Component } from "react";
import About from "./about";
import Tutorial from "./tutorial";
import AppLand from "../landing_api/Landing_page";
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
      <div style={{display: "flex", minHeight:"100vh", flexDirection: "column", justifyContent: "stretch"}}>
          <nav className="navbar navbar-expand-md"
          style={navGradient}
          > 
          {/*navbar navbar-expand-md justify-content-center m-0 navbar-dark bg-dark */}
              {/* This wraps complete navbar */}
              <div style={{padding:"0 0 0 5px", marginBottom:0, lineHeight: "50px", position:"absolute", bottom:'15px'}}>
                  <div
                  style= {this.state.home ? activeHome : inactiveHome}
                  >
                    <a  
                      style ={{color: "inherit", textDecoration:"inherit"}}
                      href="http://localhost:3000"
                      onClick={this.toggle_home}
                      >ENACTdb
                    </a>
                  </div>
              </div>

              <div className="w-100 mt-5">
                <ul className="nav navbar-nav w-100 justify-content-end">
                    <li className="nav-item">
                      <a
                        className="nav-link"
                        href="#"
                        onClick={this.toggle_about}
                        style= {this.state.about ? active : inactive}
                      >About</a>
                    </li>
                    <li className="nav-item">
                    <a
                        className="nav-link"
                        href="#"
                        onClick={this.toggle_name}
                        style= {this.state.nomenclature ? active : inactive}
                      >Nomenclature</a>
                    </li>
                    <li className="nav-item">
                    <a
                      className="nav-link"
                      href="#"
                      onClick={this.toggle_tutorial}
                      style= {this.state.tutorial ? active : inactive}
                    >Tutorial</a>
                    </li>
                </ul>
              </div>
          </nav>
        <div className="d-flex flex-grow-1 py-3">
          {/* style={{borderStyle: "solid", borderColor: "green"}} >  */}
          <div className="container p-1">
          {/* style={{borderStyle: "solid", borderColor: "blue"}}> */}
          {this.state.about ? <About/> : ' '}
          {this.state.nomenclature ? <Nomenclature/> : ' '}
          {this.state.home ? <AppLand/> : ' '}
          {/* {this.state.home ? <HolderFront/> : ' '} */}
          {this.state.tutorial ? <Tutorial/> : ' '}
          {console.log('InThere')}
          </div>
        </div>
        <div className="container p-0">
          <div className='text-center text-lg-left'
              style={footGradient}
              >
                <div className='text-center p-3 text-light' style={{textIndent: "50px"}}>
                    Server developed and maintained by shashibp-lab: "
                    <a className='text-light' href='https://shashibp-lab.github.io/' >
                      https://shashibp-lab.github.io/"
                    </a>
                  <br />
                  &copy; {new Date().getFullYear()} Copyright:{' '} Department of Biological Sciences, IISER Mohali, Sector 81, SAS
                  Nagar, 140306, Punjab, India
                </div>
          </div>
        </div>
        </div>
    );
  }
}

const boldStyle = {
  color: "white",
  fontWeight: "thick"
};
const navGradient = {
  // background:"linear-gradient(to bottom right, #2d545e, #12343b)"
  // background: "linear-gradient(to right, #52658f, #333a56",
  // above is bluish tint
  // background: "linear-gradient(to right, #00303f, #333a56",
  background: "#494e6b",
  background: "linear-gradient(to right, #494e6b, #6d7993",
  background: "linear-gradient(to right, #494e6b, #6e7376",
  background: 'linear-gradient(to right, #00303f, #333a56',
  // background: 'linear-gradient(to right, #49274a, #6e3667',
  // background: 'linear-gradient(to right, #1e392a, #015249',

  // this is cerulean to vermillion
  
  // style={{ background: 'linear-gradient(to right, #2d545e, #12343b' }} // blusih green 1st
  // style={{ background: 'linear-gradient(to right, rgba(102, 126, 234, 0.5), rgba(118, 75, 162, 0.5))' }}
  paddingBottom:"0px",
  marginBottom:"0px",
  lineHeight: "50px",
  // background: "transparent"
} 

const footGradient = {
  // background: "#192231"
  background: 'linear-gradient(to right, #00303f, #333a56',
  // background: 'linear-gradient(to right, #49274a, #6e3667', //purplish tint
  // background: 'linear-gradient(to right, #015249, #1e392a',

}
const active = {
  color: "#A9A9A9",
  fontWeight: "thick",
  fontSize: 30, 
  bottom:'0px'
};
const inactive = {
  color: "white", //carbon
  fontWeight: "thick",
  fontSize: 30, 
  bottom:'0px'
};
const activeHome = {
  color: "white",
  fontFamily: "monoton",
  fontWeight: "thick",
  fontSize: 70,
  padding:"0 0 0 5px"
  // position:'absolute', 
  // bottom:'0px'
};
const inactiveHome = {
  color: "#A9A9A9", //carbon
  fontWeight: "thick",
  fontFamily: "monoton",
  fontSize: 40,
  verticalAlign: "text-bottom",
  // position:'absolute', 
  // bottom:'0px'
};
const background = {
  //background: "#F7FSE6" // cant be distinguished from the normal white, no contrast actuially
  background: "#98878F",
  // background: "black"
  // minHeight: "600px",
  // maxHeight: "1260px",
  verticalAlign: "middle"

}
export default Navbar;
