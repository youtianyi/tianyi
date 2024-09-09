from flask import Flask,render_template,Blueprint,request,redirect,url_for,send_from_directory
from flask_script import Manager  ###pip install flask-script==2.0.5
from APP.views import blue
from APP.models import init_db
from flask_bootstrap import Bootstrap
from flask_sqlalchemy import SQLAlchemy
# create the app
app = Flask(__name__, static_url_path='/oncodrug/static', static_folder='/h/tianyi/oncodrug/static')
###注册蓝图
app.register_blueprint(blueprint=blue)
###添加boostrap
bootstrap = Bootstrap(app)
###导入初始数据库，修改配置
app.config["SQLALCHEMY_DATABASE_URI"] = "sqlite:///drug_combination.db"
app.config['SQLALCHEMY_TRACK_MODIFICATIONS']=False
init_db(app)
manager=Manager(app=app)

if __name__ == '__main__':
    app.run(host='0.0.0.0',debug=True,port=5011)
    manager.run()