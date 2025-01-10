from .factory import app, celery
from .routes import register_routes

register_routes(app)

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5002, debug=True)