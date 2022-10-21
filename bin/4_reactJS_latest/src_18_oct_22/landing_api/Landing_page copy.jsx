import React, { Component } from "react";
import Form from "../forms/Form";
import AppfetchOb from "./API_act_fetch";
import AppApi from "./API_handler";

class AppLand extends Component {
  constructor(props) {
    super(props);
    this.state = {
      post: {
        name: ""
      },
      jobs: [],
      showForm: true
    };
  }
  //have 3 state variables, showForm will render the form, post(name) and jobs [] unknown

  handleChange = e => {
    const { name, value } = e.target;
    console.log(name,value,'nameValue')
    this.setState(prevState => ({
      post: { ...prevState.post, [name]: value },
    }));
    
  };
  //When a value will be added to the user form, this change function will be called and then new value will be added to prevstate ('' in default firsttime), e.g. ...prevState.post is state value before change event is triggered or new value is added, [name]: value is, store the value
  handleSubmit = e => {
    e.preventDefault();
    this.setState(prevState => ({
      jobs: [prevState.post],
      showForm: false
    }));
  };

  render() {
    //const urlfetch = "http://172.16.3.146/nextrap/";
    //const urlfetch = "http://localhost:8001";
    //const urlfetch = "http://14.139.227.206/nextrap";
    const urlfetch = "http://localhost:8000";
    console.log("landingpage", urlfetch, this.state);
    if (this.state.showForm) {
      return (
        <div className="position-relative" align="center">
          <Form
            handleChange={this.handleChange}
            post={this.state.post}
            handleSubmit={this.handleSubmit}
          />
        </div>
      );
    } else {
      if (isNaN(this.state.jobs[0].name)) {
        console.log(
          "insidestring funct",
          this.state.jobs[0],
          isNaN(this.state.jobs[0])
        );
        return (
          <div className="post-container">
            {this.state.jobs.map((job, index) => (
              <div key={index}>
                <AppApi
                  key={index}
                  urlprefix={urlfetch}
                  urlfetch={`${urlfetch}/name/${job.name}/`}
                />
              </div>
            ))}
          </div>
        );
      } else {
        console.log("insidenumeric funct");
        return (
          <div className="post-container">
            {this.state.jobs.map((job, index) => (
              <div key={index}>
                <AppfetchOb
                  key={index}
                  urlfetch={`${urlfetch}/ncbid/${job.name}/`}
                />
                ;
              </div>
            ))}
          </div>
        );
      }
    }
  }
}
export default AppLand;
