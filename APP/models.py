from flask_sqlalchemy import SQLAlchemy
from flask import Flask



db= SQLAlchemy()
def init_db(app):
    db.init_app(app)

class DrugComdb(db.Model): ###继承model类
    ###您可以定义表中的各个列，并指定它们的名称、数据类型以及其他属性
    ###使用db.Columns定义一个简单的表结构
    Drug_combination_ID_num=db.Column(db.Integer, primary_key=True,autoincrement=True)
    Drug_combination_ID=db.Column(db.String(1000))
    Evidence_level=db.Column(db.String(1000))
    Prioritization_score=db.Column(db.Integer)
    Targeted_drug=db.Column(db.String(1000))
    Non_targeted_drug=db.Column(db.String(1000))
    Cancer_type_oncotree_Level2=db.Column(db.String(1000))
    Action_mutations=db.Column(db.String(1000))
    Response=db.Column(db.String(1000))
    Adverse_effect=db.Column(db.String(1000))
    Drug_dosage=db.Column(db.String(1000))
    Evidence_level_score=db.Column(db.Integer)
    Action_mutation_level_score=db.Column(db.Integer)
    Actionable_mutation_precision_score=db.Column(db.Integer)
    FDA_evidence_score=db.Column(db.Integer)
    Response_score=db.Column(db.Integer)
    Prioritization_score=db.Column(db.Integer)



class Match_VCF(db.Model):
    Drug_combination_ID_num=db.Column(db.Integer, primary_key=True,autoincrement=True)
    Drug_combination_ID=db.Column(db.String(1000))
    Evidence_level=db.Column(db.String(1000))
    Prioritization_score=db.Column(db.Integer)
    Targeted_drug=db.Column(db.String(1000))
    Non_targeted_drug=db.Column(db.String(1000))
    Cancer_type_oncotree_Level2=db.Column(db.String(1000))
    Action_mutations=db.Column(db.String(1000))
    Response=db.Column(db.String(1000))
    Adverse_effect=db.Column(db.String(1000))
    Drug_dosage=db.Column(db.String(1000))
    Evidence_level_score=db.Column(db.Integer)
    Action_mutation_level_score=db.Column(db.Integer)
    Actionable_mutation_precision_score=db.Column(db.Integer)
    FDA_evidence_score=db.Column(db.Integer)
    Response_score=db.Column(db.Integer)
    Prioritization_score=db.Column(db.Integer)






